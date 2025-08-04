function create_buses(zone_listing_dict::OrderedDict{String,Int64})
    # this function creates ACBus objects for each zone in the GENX model
    buses_dict = OrderedDict{String,ACBus}()
    for (zone_name, zone_number) in zone_listing_dict
        bus = ACBus(;
            available = true, #PSY5 field change
            number = zone_number, # assign zone number as bus number
            name = zone_name,  # assign string as bus name
            bustype = zone_number == 1 ? "REF" : "PV", # defining as generator bus (i.e., active power & voltage magnitude)
            angle = nothing,
            magnitude = nothing,
            voltage_limits = (min = 0.95, max = 1.05),
            base_voltage = 230.0,
            area = nothing, # we will define this next
            load_zone = nothing, #load zones are generally not used in PSY; areas are the norm
        )
        buses_dict[zone_name] = bus
    end
    return buses_dict
end

function create_areas(buses_dict::OrderedDict{String,ACBus})
    # create Area components and return dictionary
    areas_dict = OrderedDict{String,Area}()
    for (bus_name, bus) in buses_dict # loop through buses in the same order as buses_dict
        # instantiate an area object
        area_active = Area(
            name = bus_name,
            peak_active_power = 0.0,
            peak_reactive_power = 0.0,
            load_response = 0.0,
        ) 
        # add area to dictionary
        areas_dict[bus_name] = area_active
    end
    return areas_dict
end

function create_lines(existing_lines_df::DataFrame, candidate_lines_df::Union{DataFrame,Nothing}, sys::System)
    # this function creates Line objects for each transmission line in the GENX model
    lines_dict = OrderedDict{String,Line}()
    
    # Create a mapping of zone numbers (GENX) to bus names (PSY)
    GENX_zone_to_PSY_bus = OrderedDict{Int,String}()
    for bus in get_components(ACBus, sys) # loop through all buses in the system
        zone_num = get_number(bus)
        PSY_bus = get_name(bus)
        GENX_zone_to_PSY_bus[zone_num] = PSY_bus
    end
    
    # Sort the DataFrame by Network_Lines to ensure consistent order
    existing_lines_df = sort(existing_lines_df,:Network_Lines)
    
    for i in 1:nrow(existing_lines_df) # loop through all existing lines
        existing_cap = existing_lines_df[i, :Line_Max_Flow_MW]
        
        # Check to see if new transfer capacity has been added (default position is 0)
        new_cap = 0.0
        if candidate_lines_df !== nothing
            matching_lines = candidate_lines_df[candidate_lines_df.Line .== existing_lines_df[i, :Network_Lines], "New_Trans_Capacity"]
            if !isempty(matching_lines)
                new_cap = matching_lines[1]
            end
        end
        
        # Get the GENX zone numbers for start and end zones
        start_zone_num = existing_lines_df[i, :Start_Zone]
        end_zone_num = existing_lines_df[i, :End_Zone]
        
        # Get the corresponding PSY bus names
        start_bus_name = GENX_zone_to_PSY_bus[start_zone_num]
        end_bus_name = GENX_zone_to_PSY_bus[end_zone_num]
        
        # Retrieve the PSY bus objects
        start_bus = get_component(ACBus, sys, start_bus_name)
        end_bus = get_component(ACBus, sys, end_bus_name)
              
        # Create descriptive line name using bus names
        line_name = string(get_name(start_bus), "_to_", get_name(end_bus))
        
        # Create the line with proper bus connections
        # lines dont have a device_base so it gets normalized by the system base power
        # xfmr too
        line = Line(;
            name = line_name,
            available = true,
            active_power_flow = 0.0,
            reactive_power_flow = 0.0,
            arc = Arc(; from = start_bus, to = end_bus),
            r = 0.0,
            x = 0.0,
            b = (from = 0.0, to = 0.0),
            rating = round((existing_cap + new_cap) / get_base_power(sys), digits=4),
            angle_limits = (min = -0.7, max = 0.7),
        )
        # add line to lines_dict with Network_Lines as key and PSY line object as value
        lines_dict[string(existing_lines_df[i, :Network_Lines])] = line
    end
    return lines_dict
end

function create_area_interchanges(existing_lines_df::DataFrame, sys_base_power::Float64, lines_dict::OrderedDict{String,Line})
    # this function creates AreaInterchange objects for each transmission line in the GENX model
    area_interchanges_dict = OrderedDict{String,AreaInterchange}()    

    for i in 1:nrow(existing_lines_df) # loop through all existing lines
        # Define PSY Objects 
        network_line_number = string(existing_lines_df[i, :Network_Lines])
        line_psy_o = lines_dict[network_line_number]
        line_psy = get_name(line_psy_o) # string 
        line_rating = get_rating(line_psy_o)
        from_area_o = get_area(get_from(get_arc(line_psy_o)))
        to_area_o = get_area(get_to(get_arc(line_psy_o)))   

        # Define GENX objects
        reference_direction_flow = round(line_rating * existing_lines_df[i, :Profile_Forward], digits=3)
        counter_direction_flow = round(line_rating * existing_lines_df[i, :Profile_Reverse], digits=3)  
        
        # Create the AreaInterchange
        area_interchange = AreaInterchange(;
            name = line_psy, # string
            available = true, # boolean
            active_power_flow = 0.0, # float
            from_area = from_area_o, # PSY area object
            to_area = to_area_o, # PSY area object
            flow_limits = (
                from_to = round(reference_direction_flow / sys_base_power, digits=3),
                to_from = round(counter_direction_flow / sys_base_power, digits=3)
            )
        )
        
        # Add to dictionary
        area_interchanges_dict[line_psy] = area_interchange
    end
    
    return area_interchanges_dict
end

function create_transmission_interfaces(interface_df::DataFrame, lines_dict::OrderedDict{String,Line}, sys_base_power::Float64)
    # this function creates PSY TransmissionInterface objects for each simultaneous flow group in the GENX model
    # note we combine the CAISO_Export_Limit and CAISO_Import_Limit into a single interface (CAISO_Net_Import_Limit)
    interfaces_dict = OrderedDict{String,TransmissionInterface}()
    
    # Get unique flow groups
    sfg_constraints = unique(interface_df[!, "Simultaneous Flow Group"])
    
    # Convert direction strings in GENX to integers (forward=1; reverse=-1)
    interface_df.Direction = get.(Ref(Dict("forward" => 1, "reverse" => -1)), interface_df.Direction, "unknown")
    
    # Check if both CAISO_Export_Limit and CAISO_Import_Limit exist
    has_export = "CAISO_Export_Limit" in sfg_constraints
    has_import = "CAISO_Import_Limit" in sfg_constraints
    
    # Handle CAISO combined interface if both exist
    if has_export && has_import
        # Get data for both interfaces
        export_lines = interface_df[interface_df[!, "Simultaneous Flow Group"] .== "CAISO_Export_Limit", :]
        import_lines = interface_df[interface_df[!, "Simultaneous Flow Group"] .== "CAISO_Import_Limit", :]
        
        # Get limits
        export_limit = export_lines[1, "limit_MW"]
        import_limit = import_lines[1, "limit_MW"]
        
        # Create combined direction mapping
        combined_direction_mapping = Dict{String, Int}()
        
        # Add export lines direction mapping
        for row in eachrow(export_lines)
            line_number = string(row.Line_Number)
            if haskey(lines_dict, line_number)
                line_name = get_name(lines_dict[line_number])
                direction = row.Direction
                combined_direction_mapping[line_name] = direction
            else
                @warn "Line number $line_number not found in lines_dict for CAISO_Export_Limit"
            end
        end
        
        # Add import lines direction mapping
        for row in eachrow(import_lines)
            line_number = string(row.Line_Number)
            if haskey(lines_dict, line_number)
                line_name = get_name(lines_dict[line_number])
                direction = row.Direction
                combined_direction_mapping[line_name] = direction
            else
                @warn "Line number $line_number not found in lines_dict for CAISO_Import_Limit"
            end
        end
        
        # Create combined interface
        combined_interface = TransmissionInterface(
            name = "CAISO_Net_Import_Limit",
            available = true,
            active_power_flow_limits = (min = -export_limit/sys_base_power, max = import_limit/sys_base_power),
            violation_penalty = 5000.0,
            direction_mapping = combined_direction_mapping
        )
        
        # Add to dictionary
        interfaces_dict["CAISO_Net_Import_Limit"] = combined_interface
        
        # Remove both from processing list
        sfg_constraints = filter(x -> x != "CAISO_Export_Limit" && x != "CAISO_Import_Limit", sfg_constraints)
    end
    
    for interface_name in sfg_constraints # loop through all remaining simultaneous flow groups
        
        # testing with CAISO export limit
        #interface_name = "CAISO_Export_Limit"
         
        # filter df for active sfg 
        interface_lines = interface_df[interface_df[!, "Simultaneous Flow Group"] .== interface_name, :]
        
#=         println("Interface lines for $interface_name:")
        println(interface_lines)
        println("\nAvailable lines in lines_dict:")
        for (num, line) in lines_dict
            println("Line number: $num, Line name: $(get_name(line))")
        end =#
        
        # Create working dictionary to map TX line to direction (as it applies to the sfg)
        direction_mapping = Dict{String, Int}()
        for row in eachrow(interface_lines) # loop through all lines in the sfg
            # Get the PSY line name from our lines_dict using the line number
            line_number = string(row.Line_Number)
            # println("\nProcessing line number: $line_number")
            if haskey(lines_dict, line_number)
                line_name = get_name(lines_dict[line_number])
                direction = row.Direction
                direction_mapping[line_name] = direction
                # println("Added mapping: $line_name => $direction")
            else
                @warn "Line number $line_number not found in lines_dict for interface $interface_name"
            end
        end
        
        # Get the flow limit for this interface
        sfg_limit = interface_lines[1, "limit_MW"]
        
        # Create the transmission interface
        interface = TransmissionInterface(
            name = interface_name,
            available = true,
            active_power_flow_limits = (min = -9999.0/sys_base_power, max = sfg_limit/sys_base_power),
            violation_penalty = 5000.0, # making this a soft constraint for now
            direction_mapping = direction_mapping
        )
        
        # Add to dictionary
        interfaces_dict[interface_name] = interface
    end    
    return interfaces_dict
end

"""
When you define components that aren't attached to a System yet, you must define all fields related to power
in per-unit using the base_power of the component (with the exception of base_power itself, which is in MVA).
"""
function create_power_loads(demand_df::DataFrame, sys::System)
    # Create dictionary to store PowerLoad objects
    power_loads_dict = OrderedDict{String,PowerLoad}()
    
    # Create PowerLoad objects for each bus in the system
    for bus in get_components(ACBus, sys)
        # Get the bus name
        bus_name = get_name(bus)
        
        # Calculate base_power as maximum demand across ALL weather years
        # Note: column name includes "_Demand" suffix
        max_power = round(maximum(demand_df[!, Symbol(bus_name * "_Demand")]), digits=2)
        
        # Handle case where max_power is 0
        if max_power == 0.0
            @warn "Zero demand detected for bus $bus_name, setting base_power to 1.0 and per-unit attributes to 0.0"
            base_power = 1.0
            active_power = 0.0
            reactive_power = 0.0
            max_active_power = 0.0
        else
            base_power = max_power
            active_power = 0.0
            reactive_power = 0.0
            max_active_power = 1.0
        end
        
        # Create the PowerLoad object
        power_load = PowerLoad(;
            name = "Load_" * bus_name,
            available = true,
            bus = bus, # assign the bus object to the power load
            active_power = active_power,  # unitized by DEVICE base_power 
            reactive_power = reactive_power,  # unitized by DEVICE base_power 
            base_power = base_power, # MVA
            max_active_power = max_active_power,  # unitized by DEVICE base_power 
            max_reactive_power = 0.0  # No reactive power limits
        )
        
        # Add to dictionary
        power_loads_dict[bus_name] = power_load
    end
    
    return power_loads_dict
end

function create_ThermalStandard_objects(sys::System, thermal_df::DataFrame, capacity_df::DataFrame, PM_type_dict::Dict, fuel_mapping_df::DataFrame, zone_dict::OrderedDict{String,Int64})
    
    #initialize thermal standards dictionary
    thermal_standards = Dict{String, ThermalStandard}()

    for i in 1:count(!ismissing, thermal_df[:, "Resource"]) # loop through all thermal resources

        # retrieve name of thermal resource
        resource_name = thermal_df[i, :Resource]

         # Print current resource being processed
        @info "Processing thermal resource: $resource_name"

        # retrieve capacity of thermal resource
        capacity_mw = round(capacity_df[capacity_df.Resource .== resource_name, :][1, :EndCap], digits=3)

        # Skip if capacity is 0
        if capacity_mw == 0.0
            @info "Skipping thermal resource $resource_name due to zero capacity"
            continue
        end

        # retrieve min power of thermal resource (already in p.u.)
        min_power = round(thermal_df[i, :Min_Power], digits=3)

        # retrieve the zone number (GENX)
        zone_number = thermal_df[i, :Zone]
        # retrieve the bus name (PSY) by finding the key in zones_dict that matches our zone number
        bus_name = findfirst(x -> x == zone_number, zone_dict)
        if bus_name === nothing
            @error "No bus found for zone $zone_number in zones_dict"
            continue
        end
        # retrieve the bus object (PSY)
        bus_object = get_component(ACBus, sys, bus_name)

        # define heat rate curve
        if "Heat_Rate_MMBTU_per_MWh" in names(thermal_df)
            # constant avg rate HR curve
            heat_rate_curve = LinearCurve(
                round(thermal_df[i, :Heat_Rate_MMBTU_per_MWh], digits=2),
                0.0)
        else
            # constant marginal rate HR curve
            heat_rate_curve = LinearCurve(
                round(thermal_df[i, :PWFU_Heat_Rate_MMBTU_per_MWh_1], digits=2),
                round(thermal_df[i, :PWFU_Fuel_Usage_Zero_Load_MMBTU_per_h], digits=2))
        end

        # check for piece-wise linear fuel curve
        if "PWFU_Heat_Rate_MMBTU_per_MWh_2" in names(thermal_df)
            throw(ErrorException("More than one PWFU segment not supported."))
        end

        # define Variable O&M
        VOM = round(thermal_df[i, :Var_OM_Cost_per_MWh], digits=2)

        # define fixed O&M
        FOM = round(thermal_df[i, :Fixed_OM_Cost_per_MWyr], digits=2)

        # define start-up cost
        Start_Cost = round(thermal_df[i, :Start_Cost_per_MW], digits=2)

        # define fuel curve, fuel price, and VOM
        fuel_curve = FuelCurve(
            value_curve = heat_rate_curve,
            fuel_cost = 3.0, # TO-DO add in time series for varying fuel prices
            vom_cost = LinearCurve(VOM, 0.0)
        )

        # define operation cost
        Op_Cost = ThermalGenerationCost(
            variable = fuel_curve,
            fixed = FOM*capacity_mw,
            start_up = Start_Cost*capacity_mw, #TO-DO come back for fuel-related start costs
            shut_down = 0.0
        )

        # ramp limits
        ramp_up = round(thermal_df[i, :Ramp_Up_Percentage], digits=3) # GENX units %MW/hr
        ramp_down = round(thermal_df[i, :Ramp_Dn_Percentage], digits=3) # GENX units %MW/hr

        # minimum up and down time
        MUT = round(thermal_df[i, :Up_Time], digits=1) # GENX units hr
        MDT = round(thermal_df[i, :Down_Time], digits=1) # GENX units hr

        # define prime mover type using the key-value mapping from MoverTypesMapping.csv
        if !haskey(PM_type_dict, resource_name)
            @warn "No prime mover type mapping found for resource $resource_name in MoverTypesMapping.csv"
            PM_type = "OT"  # Default to Other if not found
        else
            PM_type = getproperty(PrimeMovers, Symbol(PM_type_dict[resource_name]))
        end

        # define fuel type
        fuel_GENX = thermal_df[i, :Fuel]
        matching_rows = fuel_mapping_df[fuel_mapping_df.genx_fuel_type .== fuel_GENX, :sienna_fuel_type]
        if isempty(matching_rows)
            @warn "No fuel mapping found for GENX fuel type $fuel_GENX"
            fuel_PSY = "OTHER" # assign OTHER as default
        else
            fuel_PSY = matching_rows[1]
        end

        # define thermal standard
        thermal = ThermalStandard(;
            name = resource_name,
            available = true,
            status = true,
            bus = bus_object,
            active_power = 0.0, # unitized by DEVICE base_power 
            reactive_power = 0.0, # unitized by DEVICE base_power 
            rating = 1.0, # unitized by DEVICE base_power
            active_power_limits = (min = round(min_power, digits=3), max = 1.0),
            reactive_power_limits = nothing,
            ramp_limits = (up = round(ramp_up/60, digits=3), down = round(ramp_down/60, digits=3)), # Sienna units: MW/Min
            operation_cost = Op_Cost, # TO-DO  time series for varying fuel prices and fuel-related start costs
            base_power = capacity_mw, # setting base power equal to nameplate capacity
            time_limits = (up = MUT, down = MDT), # Hours, unaffected by per-unitization
            must_run = false, # Default to false, will be updated later
            prime_mover_type = PM_type, # assign Prime Mover Type 
            fuel = fuel_PSY, #assign ThermalFuels via direct string (i.e. $FuelType - not "ThermalFuels.$FuelType")
        )

        # add thermal standard to thermal_standards dictionary  
        thermal_standards[resource_name] = thermal

    end
    return thermal_standards
end

function update_TS_must_run_status!(TS_collection::Vector{ThermalStandard})
    # Loop through all thermal standard generators
    for gen in TS_collection


        # Get the fuel type and prime mover type
        fuel_type = string(get_fuel(gen))
        pm_type = get_prime_mover_type(gen)
        
        # default value
        must_run = false # default value

        # Set must_run flag based on fuel type and prime mover type
        if fuel_type in ["WOOD_WASTE_SOLIDS", "GEOTHERMAL", "NUCLEAR"]
            must_run = true
        end
        if pm_type == PrimeMovers.OT
            must_run = true
        end
        
        # Update the must_run status
        set_must_run!(gen, must_run)
    end
end

function create_VRE_objects(sys::System, vre_df::DataFrame, capacity_df::DataFrame, PM_type_dict::Dict, zone_dict::OrderedDict{String,Int64})
    
    #initialize renewable dispatch dictionary
    RenewableDispatch_dict = Dict{String, RenewableDispatch}()

    for i in 1:count(!ismissing, vre_df[:, "Resource"]) # loop through all VRE resources

        # retrieve name of VRE resource
        resource_name = vre_df[i, :Resource]

        # Print current resource being processed
        @info "Processing VRE resource: $resource_name"

        # retrieve capacity of VRE resource
        capacity_mw = round(capacity_df[capacity_df.Resource .== resource_name, :][1, :EndCap], digits=3)

        # Skip if capacity is 0
        if capacity_mw == 0.0
            @info "Skipping VRE resource $resource_name due to zero capacity"
            continue
        end

        # retrieve the zone number (GENX)
        zone_number = vre_df[i, :Zone]
        # retrieve the bus name (PSY) by finding the key in zones_dict that matches our zone number
        bus_name = findfirst(x -> x == zone_number, zone_dict)
        if bus_name === nothing
            @error "No bus found for zone $zone_number in zones_dict"
            continue
        end
        # retrieve the bus object (PSY)
        bus_object = get_component(ACBus, sys, bus_name)

        # define Variable O&M
        VOM = round(vre_df[i, :Var_OM_Cost_per_MWh], digits=2)

        # define operation cost (simpler than thermal - just variable O&M)
        Op_Cost = RenewableGenerationCost(CostCurve(LinearCurve(VOM, 0.0)))

        # define prime mover type using the key-value mapping from MoverTypesMapping.csv
        if !haskey(PM_type_dict, resource_name)
            @warn "No prime mover type mapping found for resource $resource_name in MoverTypesMapping.csv"
            PM_type = "OT"  # Default to Other if not found
        else
            PM_type = getproperty(PrimeMovers, Symbol(PM_type_dict[resource_name]))
        end

        # define renewable dispatch
        renewable = RenewableDispatch(;
            name = resource_name,
            available = true,
            bus = bus_object,
            active_power = 0.0, # unitized by DEVICE base_power 
            reactive_power = 0.0, 
            rating = 1.0, # unitized by DEVICE base_power
            prime_mover_type = PM_type,
            reactive_power_limits = nothing, # No reactive power limits
            power_factor = 1.0, # Unity power factor
            operation_cost = Op_Cost,
            base_power = round(capacity_mw, digits=3), # setting base power to nameplate capacity; DONT FORGET THE COMMA
        )

        # add renewable dispatch to vre_standards dictionary  
        RenewableDispatch_dict[resource_name] = renewable

    end
    return RenewableDispatch_dict
end

function create_btm_objects(sys::System, vre_df::DataFrame, capacity_df::DataFrame, PM_type_dict::Dict, zone_dict::OrderedDict{String,Int64})
    
    #initialize renewable dispatch dictionary
    RenewableNonDispatch_dict = Dict{String, RenewableNonDispatch}()

    for i in 1:count(!ismissing, btm_df[:, "Resource"]) # loop through all VRE resources

        # retrieve name of VRE resource
        resource_name = btm_df[i, :Resource]

        # Print current resource being processed
        @info "Processing BTM resource: $resource_name"

        # retrieve capacity of BTM resource
        capacity_mw = round(capacity_df[capacity_df.Resource .== resource_name, :][1, :EndCap], digits=3)

        # Skip if capacity is 0
        if capacity_mw == 0.0
            @info "Skipping BTM resource $resource_name due to zero capacity"
            continue
        end

        # retrieve the zone number (GENX)
        zone_number = btm_df[i, :Zone]
        # retrieve the bus name (PSY) by finding the key in zones_dict that matches our zone number
        bus_name = findfirst(x -> x == zone_number, zone_dict)
        if bus_name === nothing
            @error "No bus found for zone $zone_number in zones_dict"
            continue
        end
        # retrieve the bus object (PSY)
        bus_object = get_component(ACBus, sys, bus_name)

        # define Variable O&M
        VOM = round(vre_df[i, :Var_OM_Cost_per_MWh], digits=2)

        # define fixed O&M
        FOM = round(btm_df[i, :Fixed_OM_Cost_per_MWyr], digits=2)

        # define operation cost (simpler than thermal - just variable O&M)
        Op_Cost = RenewableGenerationCost(CostCurve(LinearCurve(VOM, 0.0)))

        # define prime mover type using the key-value mapping from MoverTypesMapping.csv
        if !haskey(PM_type_dict, resource_name)
            @warn "No prime mover type mapping found for resource $resource_name in MoverTypesMapping.csv"
            PM_type = PrimeMovers.OT  # Default to Other if not found
        else
            PM_type = getproperty(PrimeMovers, Symbol(PM_type_dict[resource_name]))
        end

        # define renewable dispatch
        renewable_ND = RenewableNonDispatch(;
            name = resource_name,
            available = true,
            bus = bus_object,
            active_power = 0.0, # unitized by DEVICE base_power 
            reactive_power = 0.0, # unitized by DEVICE base_power 
            rating = 1.0, # unitized by DEVICE base_power
            prime_mover_type = PM_type,
            power_factor = 1.0, # Unity power factor
            # operation_cost = Op_Cost, # not available for RenewableNonDispatch
            base_power = capacity_mw, # setting base power to nameplate capacity
        )

        # add renewable dispatch to vre_standards dictionary  
        RenewableNonDispatch_dict[resource_name] = renewable_ND

    end
    return RenewableNonDispatch_dict
end

function create_Hydro_objects(sys::System, hydro_df::DataFrame, capacity_df::DataFrame, PM_type_dict::Dict, zone_dict::OrderedDict{String,Int64})
    
    #initialize hydro dispatch dictionary
    HydroDispatch_dict = Dict{String, HydroDispatch}()

    for i in 1:count(!ismissing, hydro_df[:, "Resource"]) # loop through all hydro resources

        # retrieve name of hydro resource
        resource_name = hydro_df[i, :Resource]

        # Print current resource being processed
        @info "Processing hydro resource: $resource_name"

        # retrieve capacity of hydro resource
        capacity_mw = round(capacity_df[capacity_df.Resource .== resource_name, :][1, :EndCap], digits=3)

        # Skip if capacity is 0
        if capacity_mw == 0.0
            @info "Skipping hydro resource $resource_name due to zero capacity"
            continue
        end

        # retrieve min power of hydro resource (already in p.u.)
        min_power = round(hydro_df[i, :Min_Power], digits=3)

        # retrieve the zone number (GENX)
        zone_number = hydro_df[i, :Zone]
        # retrieve the bus name (PSY) by finding the key in zones_dict that matches our zone number
        bus_name = findfirst(x -> x == zone_number, zone_dict)
        if bus_name === nothing
            @error "No bus found for zone $zone_number in zones_dict"
            continue
        end
        # retrieve the bus object (PSY)
        bus_object = get_component(ACBus, sys, bus_name)

        # define Variable O&M
        VOM = round(hydro_df[i, :Var_OM_Cost_per_MWh], digits=2)

        # define fixed O&M
        FOM = round(hydro_df[i, :Fixed_OM_Cost_per_MWyr], digits=2)

        # define operation cost (simpler than thermal - just variable O&M)
        Op_Cost = HydroGenerationCost(
            variable = CostCurve(LinearCurve(VOM, 0.0)), # ProductionVariableCostCurve
            fixed = FOM*capacity_mw # float
        )

        # ramp limits
        ramp_up = hydro_df[i, :Ramp_Up_Percentage] # GENX units %MW/hr
        ramp_down = hydro_df[i, :Ramp_Dn_Percentage] # GENX units %MW/hr
        ramp_up_sienna = round(ramp_up/60, digits=3) # Sienna units: MW/Min
        ramp_down_sienna = round(ramp_down/60, digits=3) # Sienna units: MW/Min

        # define prime mover type using the key-value mapping from MoverTypesMapping.csv
        if !haskey(PM_type_dict, resource_name)
            @warn "No prime mover type mapping found for resource $resource_name in MoverTypesMapping.csv"
            PM_type = PrimeMovers.OT  # Default to Other if not found
        else
            PM_type = getproperty(PrimeMovers, Symbol(PM_type_dict[resource_name]))
        end

        # define hydro dispatch (i.e., run-of-river)
        hydro = HydroDispatch(;
            name = resource_name,
            available = true,
            bus = bus_object,
            active_power = 0.0, # unitized by DEVICE base_power 
            reactive_power = 0.0, # unitized by DEVICE base_power 
            rating = 1.0, # unitized by DEVICE base_power
            prime_mover_type = PM_type,
            active_power_limits = (min = min_power, max = 1.0),
            reactive_power_limits = nothing, # No reactive power limits
            ramp_limits = (up = ramp_up_sienna, down = ramp_down_sienna), # Sienna units: MW/Min
            time_limits = nothing, # units: Hours -> not defined in GenX
            operation_cost = Op_Cost,
            base_power = capacity_mw, # setting base power to nameplate capacity
        )

        # add hydro dispatch to HydroDispatch_dict  
        HydroDispatch_dict[resource_name] = hydro

    end
    return HydroDispatch_dict
end

function create_storage_objects(sys::System, storage_df::DataFrame, capacity_df::DataFrame, PM_type_dict::Dict, storage_type_dict::Dict, zone_dict::OrderedDict{String,Int64})

    #initialize storage dictionary
    Storage_dict = Dict{String, EnergyReservoirStorage}()

    for i in 1:count(!ismissing, storage_df[:, "Resource"]) # loop through all storage resources

        # retrieve name of storage resource
        resource_name = storage_df[i, :Resource]

        # Print current resource being processed
        @info "Processing storage resource: $resource_name"

        # Skip if it's a pumped hydro unit
        if haskey(PM_type_dict, resource_name) && PM_type_dict[resource_name] == "PS"
            @info "Skipping pumped hydro resource $resource_name - will be handled separately"
            continue
        end

        # retrieve capacity of storage resource
        power_capacity_mw = round(capacity_df[capacity_df.Resource .== resource_name, :][1, :EndCap], digits=3)
        energy_capacity_mwh = round(capacity_df[capacity_df.Resource .== resource_name, :][1, :EndEnergyCap], digits=3) 

        # Skip if either power or energy capacity is 0
        if power_capacity_mw == 0.0 || energy_capacity_mwh == 0.0
            @info "Skipping storage resource $resource_name due to zero capacity (power: $power_capacity_mw MW, energy: $energy_capacity_mwh MWh)"
            continue
        end

        # retrieve the zone number (GENX)
        zone_number = storage_df[i, :Zone]
        # retrieve the bus name (PSY) by finding the key in zones_dict that matches our zone number
        bus_name = findfirst(x -> x == zone_number, zone_dict)
        if bus_name === nothing
            @error "No bus found for zone $zone_number in zones_dict"
            continue
        end
        # retrieve the bus object (PSY)
        bus_object = get_component(ACBus, sys, bus_name)

        # define Variable O&M costs
        charge_VOM = round(storage_df[i, :Var_OM_Cost_per_MWh], digits=2)
        discharge_VOM = round(storage_df[i, :Var_OM_Cost_per_MWh], digits=2)

        # Fixed O&M
        FOM_MW = round(storage_df[i, :Fixed_OM_Cost_per_MWyr], digits=2)
        FOM_MWh = round(storage_df[i, :Fixed_OM_Cost_per_MWhyr], digits=2)

        # define operation cost
        Op_Cost = StorageCost(
            charge_variable_cost = CostCurve(LinearCurve(charge_VOM, 0.0)),
            discharge_variable_cost = CostCurve(LinearCurve(discharge_VOM, 0.0)),
            fixed = FOM_MW*power_capacity_mw + FOM_MWh*energy_capacity_mwh,
            start_up = 0.0,
            shut_down = 0.0,
            energy_shortage_cost = 0.0, # Cost incurred by the model for being short of the energy target
            energy_surplus_cost = 0.0, # Cost incurred by the model for surplus energy stored
        )

        # define prime mover type using the key-value mapping from MoverTypesMapping.csv
        if !haskey(PM_type_dict, resource_name)
            @warn "No prime mover type mapping found for resource $resource_name in MoverTypesMapping.csv"
            PM_type = PrimeMovers.OT  # Default to Other if not found
        else
            PM_type = getproperty(PrimeMovers, Symbol(PM_type_dict[resource_name]))
        end

        # define storage technology type
        if !haskey(storage_type_dict, resource_name)
            @warn "No storage technology type mapping found for resource $resource_name in StorageMapping.csv"
            ST_type = StorageTech.BAT  # Default to Battery if not found
        else
            ST_type = getproperty(StorageTech, Symbol(storage_type_dict[resource_name]))
        end

        # define storage efficiency
        charge_efficiency = round(storage_df[i, :Eff_Up], digits=3)
        discharge_efficiency = round(storage_df[i, :Eff_Down], digits=3)

        # define storage device
        storage = EnergyReservoirStorage(;
            name = resource_name,
            available = true,
            bus = bus_object,
            prime_mover_type = PM_type,
            storage_technology_type = ST_type,
            storage_capacity = power_capacity_mw == 0 ? 0.0 : round(energy_capacity_mwh/power_capacity_mw, digits=3),# unitized by device base power 
            storage_level_limits = (min = 0.0, max = 1.0), # limits on SOC range 
            initial_storage_capacity_level = 0.50, # initial SOC level 
            rating = 1.0, # max output power rating; unitized by DEVICE base_power
            active_power = 0.0, # initial active power output 
            input_active_power_limits = (min = 0.0, max = 1.0),
            output_active_power_limits = (min = 0.0, max = 1.0),
            efficiency = (in = charge_efficiency, out = discharge_efficiency),
            reactive_power = 0.0, # unitized by DEVICE base_power 
            reactive_power_limits = nothing, # No reactive power limits
            base_power = power_capacity_mw, # setting base power to nameplate capacity
            operation_cost = Op_Cost,
            conversion_factor = 1.0, # Conversion factor of storage_capacity to MWh, if different than 1.0.
            storage_target = 0.0, #  Storage target at the end of simulation as ratio of storage capacity
            cycle_limits = 365, # Storage Maximum number of cycles per year
        )

        # add storage device to Storage_dict  
        Storage_dict[resource_name] = storage

    end
    return Storage_dict
end

function create_PHS_storage_objects(sys::System, storage_df::DataFrame, capacity_df::DataFrame, PM_type_dict::Dict, storage_type_dict::Dict, zone_dict::OrderedDict{String,Int64})

    #initialize storage dictionary
    PHS_Storage_dict = Dict{String, EnergyReservoirStorage}()

    for i in 1:count(!ismissing, storage_df[:, "Resource"]) # loop through all storage resources (Batteries and PHS)

        # retrieve name of storage resource
        resource_name = storage_df[i, :Resource]

        # Skip if it's a pumped hydro unit
        if haskey(PM_type_dict, resource_name) && PM_type_dict[resource_name] == "PS"
            # Print current resource being processed
            @info "Processing PHS storage resource: $resource_name"
        else
            continue # terminate current iteration and proceed to next resource
        end

        # retrieve capacity of storage resource
        power_capacity_mw = round(capacity_df[capacity_df.Resource .== resource_name, :][1, :EndCap], digits=3)
        energy_capacity_mwh = round(capacity_df[capacity_df.Resource .== resource_name, :][1, :EndEnergyCap], digits=3) 

        # Skip if either power or energy capacity is 0
        if power_capacity_mw == 0.0 || energy_capacity_mwh == 0.0
            @info "Skipping storage resource $resource_name due to zero capacity (power: $power_capacity_mw MW, energy: $energy_capacity_mwh MWh)"
            continue
        end

        # retrieve the zone number (GENX)
        zone_number = storage_df[i, :Zone]
        # retrieve the bus name (PSY) by finding the key in zones_dict that matches our zone number
        bus_name = findfirst(x -> x == zone_number, zone_dict)
        if bus_name === nothing
            @error "No bus found for zone $zone_number in zones_dict"
            continue
        end
        # retrieve the bus object (PSY)
        bus_object = get_component(ACBus, sys, bus_name)

        # define Variable O&M costs
        charge_VOM = round(storage_df[i, :Var_OM_Cost_per_MWh], digits=2)
        discharge_VOM = round(storage_df[i, :Var_OM_Cost_per_MWh], digits=2)

        # Fixed O&M
        FOM_MW = round(storage_df[i, :Fixed_OM_Cost_per_MWyr], digits=2)
        FOM_MWh = round(storage_df[i, :Fixed_OM_Cost_per_MWhyr], digits=2)

        # define operation cost
        Op_Cost = StorageCost(
            charge_variable_cost = CostCurve(LinearCurve(charge_VOM, 0.0)),
            discharge_variable_cost = CostCurve(LinearCurve(discharge_VOM, 0.0)),
            fixed = FOM_MW*power_capacity_mw + FOM_MWh*energy_capacity_mwh,
            start_up = 0.0,
            shut_down = 0.0,
            energy_shortage_cost = 0.0, # Cost incurred by the model for being short of the energy target
            energy_surplus_cost = 0.0, # Cost incurred by the model for surplus energy stored
        )

        # define prime mover type using the key-value mapping from MoverTypesMapping.csv
        if !haskey(PM_type_dict, resource_name)
            @warn "No prime mover type mapping found for resource $resource_name in MoverTypesMapping.csv"
            PM_type = PrimeMovers.OT  # Default to Other if not found
        else
            PM_type = getproperty(PrimeMovers, Symbol(PM_type_dict[resource_name]))
        end

        # define storage technology type
        if !haskey(storage_type_dict, resource_name)
            @warn "No storage technology type mapping found for resource $resource_name in StorageMapping.csv"
            ST_type = StorageTech.BAT  # Default to Battery if not found
        else
            ST_type = getproperty(StorageTech, Symbol(storage_type_dict[resource_name]))
        end

        # define storage efficiency
        charge_efficiency = round(storage_df[i, :Eff_Up], digits=3)
        discharge_efficiency = round(storage_df[i, :Eff_Down], digits=3)

        # define PHS storage device (using battery storage as our template)
        PHS_storage = EnergyReservoirStorage(;
            name = resource_name,
            available = true,
            bus = bus_object,
            prime_mover_type = PM_type,
            storage_technology_type = ST_type,
            storage_capacity = power_capacity_mw == 0 ? 0.0 : round(energy_capacity_mwh/power_capacity_mw, digits=3),# unitized by device base power 
            storage_level_limits = (min = 0.0, max = 1.0), # limits on SOC range 
            initial_storage_capacity_level = 0.50, # initial SOC level 
            rating = 1.0, # max output power rating; unitized by DEVICE base_power
            active_power = 0.0, # initial active power output 
            input_active_power_limits = (min = 0.0, max = 1.0),
            output_active_power_limits = (min = 0.0, max = 1.0),
            efficiency = (in = charge_efficiency, out = discharge_efficiency),
            reactive_power = 0.0, # unitized by DEVICE base_power 
            reactive_power_limits = nothing, # No reactive power limits
            base_power = power_capacity_mw, # setting base power to nameplate capacity
            operation_cost = Op_Cost,
            conversion_factor = 1.0, # Conversion factor of storage_capacity to MWh, if different than 1.0.
            storage_target = 0.0, #  Storage target at the end of simulation as ratio of storage capacity
            cycle_limits = 365, # Storage Maximum number of cycles per year
        )

        # add storage device to Storage_dict  
        PHS_Storage_dict[resource_name] = PHS_storage

    end
    return PHS_Storage_dict
end

#= function create_PHS_objects(sys::System, storage_df::DataFrame, capacity_df::DataFrame, PM_type_dict::Dict, zone_dict::OrderedDict{String,Int64})
    
    #initialize pumped hydro dictionary
    PumpedHydro_dict = Dict{String, HydroPumpedStorage}()

    for i in 1:count(!ismissing, storage_df[:, "Resource"]) #note: PHS parameters are in Storage.csv
        # retrieve name of storage resource
        resource_name = storage_df[i, :Resource]

        # Skip if it's not a pumped hydro unit
        if !haskey(PM_type_dict, resource_name) || PM_type_dict[resource_name] != "PS"
            continue
        end

        # Print current resource being processed
        @info "Processing pumped hydro resource: $resource_name"

        # retrieve capacity of pumped hydro resource
        power_capacity_mw = round(capacity_df[capacity_df.Resource .== resource_name, :][1, :EndCap], digits=3)
        energy_capacity_mwh = round(capacity_df[capacity_df.Resource .== resource_name, :][1, :EndEnergyCap], digits=3)

        # Skip if either power or energy capacity is 0
        if power_capacity_mw == 0.0 || energy_capacity_mwh == 0.0
            @info "Skipping pumped hydro resource $resource_name due to zero capacity (power: $power_capacity_mw MW, energy: $energy_capacity_mwh MWh)"
            continue
        end

        # retrieve the zone number (GENX)
        zone_number = storage_df[i, :Zone]
        # retrieve the bus name (PSY) by finding the key in zones_dict that matches our zone number
        bus_name = findfirst(x -> x == zone_number, zone_dict)
        if bus_name === nothing
            @error "No bus found for zone $zone_number in zones_dict"
            continue
        end
        # retrieve the bus object (PSY)
        bus_object = get_component(ACBus, sys, bus_name)

        # define Variable O&M costs
        charge_VOM = round(storage_df[i, :Var_OM_Cost_per_MWh], digits=2)
        discharge_VOM = round(storage_df[i, :Var_OM_Cost_per_MWh], digits=2)

        # Fixed O&M
        FOM_MW = round(storage_df[i, :Fixed_OM_Cost_per_MWyr], digits=2)
        FOM_MWh = round(storage_df[i, :Fixed_OM_Cost_per_MWhyr], digits=2)

        # define operation cost
        Op_Cost = StorageCost(
            charge_variable_cost = CostCurve(LinearCurve(charge_VOM, 0.0)),
            discharge_variable_cost = CostCurve(LinearCurve(discharge_VOM, 0.0)),
            fixed = FOM_MW*power_capacity_mw + FOM_MWh*energy_capacity_mwh,
            start_up = 0.0,
            shut_down = 0.0,
            energy_shortage_cost = 0.0,
            energy_surplus_cost = 0.0,
        )

        # define storage efficiency
        charge_efficiency = round(storage_df[i, :Eff_Up], digits=3) #note: currently  applied for both charge and discharge
        # discharge_efficiency = round(storage_df[i, :Eff_Down], digits=3) #not used in PSY's PHS formulation 
        PHS_efficiency = round(sqrt(charge_efficiency), digits=3) # half this to get proper roundtrip efficiency

        # calculate duration in hours
        duration_hours = round(energy_capacity_mwh/(power_capacity_mw*PHS_efficiency), digits=3)
        # calculate initial storage (50% of reservoir)
        initial_storage_hours = round(duration_hours * 0.5, digits=3)

        # define pumped hydro storage device
        pumped_hydro = HydroPumpedStorage(;
            name = resource_name, # string
            available = true, # Boolean
            bus = bus_object, # PSY bus object 
            active_power = 0.0, # initial active power output
            reactive_power = 0.0, # initial reactive power output
            rating = 1.0, # max output power rating; unitized by DEVICE base_power
            base_power = power_capacity_mw, # setting base power equal to nameplate capacity
            prime_mover_type = PrimeMovers.PS, # Pumped Storage (hard coded)
            active_power_limits = (min = 0.0, max = 1.0), # generation mode limits
            reactive_power_limits = nothing, # No reactive power limits
            ramp_limits = nothing, # No ramp limits
            time_limits = nothing, # No time limits
            rating_pump = 1.0, # max input power rating; unitized by DEVICE base_power
            active_power_limits_pump = (min = 0.0, max = 1.0), # pumping mode limits
            reactive_power_limits_pump = nothing, # No reactive power limits in pump mode
            ramp_limits_pump = nothing, # No ramp limits (not defined for pumped hydro)
            time_limits_pump = nothing, # No time limits (not defined for pumped hydro)
            storage_capacity = (up = duration_hours, down = duration_hours*2), # setting lower reservoir equal to 2x the upper reservoir
            inflow = 0.0, # mandatory object; we will assign later
            outflow = 0.0, # mandatory object; we will assign later
            initial_storage = (up = initial_storage_hours, down = initial_storage_hours*2), # initial storage level; units: hours
            storage_target = (up = 0.5, down = 0.5), # Storage target of UPPER res @ end of simulation (ratio of storage capacity)
            operation_cost = Op_Cost,
            pump_efficiency = PHS_efficiency, # roundtrip efficiency
            conversion_factor = 1.0, # Conversion factor from flow to energy p.u.-hr
            status = PSY.PumpHydroStatusModule.PumpHydroStatus.OFF, # initial status; default: "PumpHydroStatus.OFF"
            time_at_status = 10.0, # initial time at status
            services = Device[], # no services
            dynamic_injector = nothing, # no dynamic injector
        )

        # add pumped hydro device to PumpedHydro_dict  
        PumpedHydro_dict[resource_name] = pumped_hydro
    end

    return PumpedHydro_dict
end =#




function create_CAISO_reg_reserve_services()
    # Initialize reserve dictionary
    Reserve_dict = Dict{String, Reserve}()
    
    # Create the VariableReserve object for up reserves
    reserve_up = VariableReserve{ReserveUp}(;
        name = "CAISO_reg_up", # string
        available = true, # Boolean representing if the reserve is active or not
        time_frame = 10.0,  # saturation timeframe to provide service (min)
        requirement = 1.0,  # scaled by timeseries (which we attach later)
        sustained_time = 3600.0,  # time (secs) reserve contribution must be sustained
        max_output_fraction = 1.0, #  the max fraction of each device's output that can be assigned
        max_participation_factor = 1.0, # the max portion [0, 1.0] of the reserve that can be contributed per device
        deployed_fraction = 0.0, # Fraction of service procurement that is assumed to be actually deployed
    )
    
    # Create the VariableReserve object for down reserves
    reserve_down = VariableReserve{ReserveDown}(;
        name = "CAISO_reg_down", # string
        available = true, # Boolean representing if the reserve is active or not
        time_frame = 10.0,  # saturation timeframe to provide service (min)
        requirement = 1.0,  # scaled by the requirement timeseries (which we attach to object later)
        sustained_time = 3600.0,  # time (secs) reserve contribution must be sustained
        max_output_fraction = 1.0, #  the max fraction of each device's output that can be assigned
        max_participation_factor = 1.0, # the max portion [0, 1.0] of the reserve that can be contributed per device
        deployed_fraction = 0.0, # Fraction of service procurement that is assumed to be actually deployed
    )
    
    # Add to dictionary
    Reserve_dict["CAISO_reg_up"] = reserve_up
    Reserve_dict["CAISO_reg_down"] = reserve_down
    
    return Reserve_dict
end

function create_CAISO_reg_reserve_units(
    thermal_generators::Vector{ThermalStandard},
    renew_d_generators::Vector{RenewableDispatch},
    hydro_generators::Vector{HydroDispatch},
    storage_units::Vector{EnergyReservoirStorage}, 
    #pumped_hydro_units::Vector{HydroPumpedStorage},
    pumped_hydro_units::Vector{EnergyReservoirStorage}, # temporary fix for PHS
    thermal_df::DataFrame,
    vre_df::DataFrame,
    hydro_df::DataFrame,
    storage_df::DataFrame
)

    # right now only assigning reserves to theral and batteries (no pumped hydro)

    # Initialize dictionary to store eligible resources
    eligible_reg_units_dict = Dict{String, Vector{Device}}()
    
    # Define utility areas we want to include
    utility_areas = ["PGE", "SCE", "SDGE"]
    
    # Initialize vectors for both up and down reserves
    eligible_reg_units_dict["CAISO_reg_up"] = Vector{Device}()
    eligible_reg_units_dict["CAISO_reg_down"] = Vector{Device}()
    
    # Check thermal generators
    for gen in thermal_generators
        bus_name = get_name(get_bus(gen))
        if bus_name in utility_areas
            # Get the resource name and check Reg_Max in the original dataframe
            resource_name = get_name(gen)
            thermal_row = thermal_df[thermal_df.Resource .== resource_name, :]
            if !isempty(thermal_row) && thermal_row[1, :Reg_Max] > 0
                # Add to both up and down reserve lists
                push!(eligible_reg_units_dict["CAISO_reg_up"], gen)
                push!(eligible_reg_units_dict["CAISO_reg_down"], gen)
            end
        end
    end
    
#=     # Check renewable dispatch generators
    for gen in renew_d_generators
        bus_name = get_name(get_bus(gen))
        if bus_name in utility_areas
            # Get the resource name and check Reg_Max in the original dataframe
            resource_name = get_name(gen)
            vre_row = vre_df[vre_df.Resource .== resource_name, :]
            if !isempty(vre_row) && vre_row[1, :Reg_Max] > 0
                # Add to both up and down reserve lists
                push!(eligible_reg_units_dict["CAISO_reg_up"], gen)
                push!(eligible_reg_units_dict["CAISO_reg_down"], gen)
            end
        end
    end =#
    
#=     # Check hydro generators
    for gen in hydro_generators
        bus_name = get_name(get_bus(gen))
        if bus_name in utility_areas
            # Get the resource name and check Reg_Max in the original dataframe
            resource_name = get_name(gen)
            hydro_row = hydro_df[hydro_df.Resource .== resource_name, :]
            if !isempty(hydro_row) && hydro_row[1, :Reg_Max] > 0
                # Add to both up and down reserve lists
                push!(eligible_reg_units_dict["CAISO_reg_up"], gen)
                push!(eligible_reg_units_dict["CAISO_reg_down"], gen)
            end
        end
    end =#
    
    # Check storage resources
    for storage in storage_units
        bus_name = get_name(get_bus(storage))
        if bus_name in utility_areas
            # Get the resource name and check Reg_Max in the original dataframe
            resource_name = get_name(storage)
            storage_row = storage_df[storage_df.Resource .== resource_name, :]
            if !isempty(storage_row) && storage_row[1, :Reg_Max] > 0
                # Add to both up and down reserve lists
                push!(eligible_reg_units_dict["CAISO_reg_up"], storage)
                push!(eligible_reg_units_dict["CAISO_reg_down"], storage)
            end
        end
    end

#=     # Check pumped hydro resources  
    for PHS in pumped_hydro_units
        bus_name = get_name(get_bus(PHS))
        if bus_name in utility_areas
            # Get the resource name and check Reg_Max in the original dataframe
            resource_name = get_name(PHS)
            PHS_row = storage_df[storage_df.Resource .== resource_name, :]
            if !isempty(PHS_row) && PHS_row[1, :Reg_Max] > 0
                # Add to both up and down reserve lists
                push!(eligible_reg_units_dict["CAISO_reg_up"], PHS)
                push!(eligible_reg_units_dict["CAISO_reg_down"], PHS)
            end
        end
    end =#
    
    return eligible_reg_units_dict
end

