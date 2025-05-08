function create_demand_PSY_timeseries(demand_ts::DataFrame, power_loads_dict::OrderedDict{String,PowerLoad})
    # Initialize container for storing timeseries as nested dictionary
    ts_container = Dict{String, Dict{String, SingleTimeSeries}}()
       
    # Define fixed calendar year timestamps (2035)
    tstamps = collect(range(DateTime("2035-01-01T00:00:00"), DateTime("2035-12-31T23:00:00"), step = Dates.Hour(1)))
    
    # Get unique years from the DateTime column and sort them
    years = sort(unique(Dates.year.(demand_ts.DateTime)))
    
    # Loop through each power load
    for (bus_name, power_load) in power_loads_dict
        # Initialize inner dictionary for this bus
        ts_container[bus_name] = Dict{String, SingleTimeSeries}()
        
        # Get the base power for normalization
        base_power = get_base_power(power_load)
        
        # Loop through each year in sorted order
        for year in years
            # Filter data for this year and sort by DateTime
            year_data = sort(demand_ts[Dates.year.(demand_ts.DateTime) .== year, :], :DateTime)
            
            # Get the demand column name for this bus
            demand_col = Symbol(bus_name * "_Demand")
            
            # Verify hourly resolution, accounting for leap day skips
            for i in 2:length(year_data.DateTime)
                time_diff = year_data.DateTime[i] - year_data.DateTime[i-1]
                is_leap_day_skip = (Dates.month(year_data.DateTime[i-1]) == 2 && Dates.day(year_data.DateTime[i-1]) == 28 &&
                                  Dates.month(year_data.DateTime[i]) == 3 && Dates.day(year_data.DateTime[i]) == 1)
                
                if time_diff != Dates.Hour(1) && !is_leap_day_skip
                    @warn "Non-hourly resolution detected in year $year for $bus_name at index $i"
                end
            end
            
            # Normalize the data by base_power
            normalized_data = year_data[!, demand_col] ./ base_power
            
            # Create the timeseries name with year suffix
            ts_name = "max_active_power_$year"
            
            # Create the SingleTimeSeries object using fixed 2035 timestamps
            ts = SingleTimeSeries(;
                name = ts_name,
                data = TimeArray(tstamps, normalized_data),
                scaling_factor_multiplier = get_max_active_power,
            )
            
            # Add to container using nested dictionary structure
            ts_container[bus_name][string(year)] = ts
        end
    end
    
    return ts_container
end

function create_Renew_D_PSY_timeseries(gen_variability_ts::DataFrame, Renew_D_collection::Vector{RenewableDispatch})
    # Initialize container for storing timeseries as nested dictionary
    ts_container = Dict{String, Dict{String, SingleTimeSeries}}()
        
    # Define fixed calendar year timestamps (2035)
    tstamps = collect(range(DateTime("2035-01-01T00:00:00"), DateTime("2035-12-31T23:00:00"), step = Dates.Hour(1)))
    
    # Get unique years from the DateTime column and sort them
    years = sort(unique(Dates.year.(gen_variability_ts.DateTime)))
    
    # Loop through each renewable dispatch generator
    #for renewable in Renew_D_collection[1:1] # testing with first element
    for renewable in Renew_D_collection
        # Get the resource name and base power
        resource_name = get_name(renewable)
        base_power = 1 # profiles are already normalized to base power
        
        # Initialize inner dictionary for this resource
        ts_container[resource_name] = Dict{String, SingleTimeSeries}()
        
        # Loop through each year in sorted order
        for year in years
            # Filter data for this year and sort by DateTime
            year_data = sort(gen_variability_ts[Dates.year.(gen_variability_ts.DateTime) .== year, :], :DateTime)
            
            # Get the availability column for this generator
            if !hasproperty(year_data, Symbol(resource_name))
                @warn "No availability data found for generator $resource_name in year $year"
                continue
            end
            
            # Verify hourly resolution, accounting for leap day skips
            for i in 2:length(year_data.DateTime)
                time_diff = year_data.DateTime[i] - year_data.DateTime[i-1]
                is_leap_day_skip = (Dates.month(year_data.DateTime[i-1]) == 2 && Dates.day(year_data.DateTime[i-1]) == 28 &&
                                  Dates.month(year_data.DateTime[i]) == 3 && Dates.day(year_data.DateTime[i]) == 1)
                
                if time_diff != Dates.Hour(1) && !is_leap_day_skip
                    @warn "Non-hourly resolution detected in year $year for $resource_name at index $i"
                end
            end
            
            # Normalize the data by base_power
            normalized_data = year_data[!, Symbol(resource_name)] ./ base_power
            
            # Create the timeseries name with year suffix
            ts_name = "max_active_power_$year"
            
            # Create the SingleTimeSeries object using fixed 2035 timestamps
            ts = SingleTimeSeries(;
                name = ts_name,
                data = TimeArray(tstamps, normalized_data),
                scaling_factor_multiplier = get_max_active_power,
            )
            
            # Add to container using nested dictionary structure
            ts_container[resource_name][string(year)] = ts
        end
    end
    
    return ts_container
end

function create_Renew_ND_PSY_timeseries(gen_variability_ts::DataFrame, Renew_ND_collection::Vector{RenewableNonDispatch})
    # Initialize container for storing timeseries as nested dictionary
    ts_container = Dict{String, Dict{String, SingleTimeSeries}}()
        
    # Define fixed calendar year timestamps (2035)
    tstamps = collect(range(DateTime("2035-01-01T00:00:00"), DateTime("2035-12-31T23:00:00"), step = Dates.Hour(1)))
    
    # Get unique years from the DateTime column and sort them
    years = sort(unique(Dates.year.(gen_variability_ts.DateTime)))
    
    # Loop through each renewable dispatch generator
    #for renewable in Renew_ND_collection[1:1] # testing with first element
    for renewable in Renew_ND_collection
        # Get the resource name and base power
        resource_name = get_name(renewable)
        base_power = 1 # profiles are already normalized to base power
        
        # Initialize inner dictionary for this resource
        ts_container[resource_name] = Dict{String, SingleTimeSeries}()
        
        # Loop through each year in sorted order
        for year in years
            # Filter data for this year and sort by DateTime
            year_data = sort(gen_variability_ts[Dates.year.(gen_variability_ts.DateTime) .== year, :], :DateTime)
            
            # Get the availability column for this generator
            if !hasproperty(year_data, Symbol(resource_name))
                @warn "No availability data found for generator $resource_name in year $year"
                continue
            end
            
            # Verify hourly resolution, accounting for leap day skips
            for i in 2:length(year_data.DateTime)
                time_diff = year_data.DateTime[i] - year_data.DateTime[i-1]
                is_leap_day_skip = (Dates.month(year_data.DateTime[i-1]) == 2 && Dates.day(year_data.DateTime[i-1]) == 28 &&
                                  Dates.month(year_data.DateTime[i]) == 3 && Dates.day(year_data.DateTime[i]) == 1)
                
                if time_diff != Dates.Hour(1) && !is_leap_day_skip
                    @warn "Non-hourly resolution detected in year $year for $resource_name at index $i"
                end
            end
            
            # Normalize the data by base_power
            normalized_data = year_data[!, Symbol(resource_name)] ./ base_power
            
            # Create the timeseries name with year suffix
            ts_name = "max_active_power_$year"
            
            # Create the SingleTimeSeries object using fixed 2035 timestamps
            ts = SingleTimeSeries(;
                name = ts_name,
                data = TimeArray(tstamps, normalized_data),
                scaling_factor_multiplier = get_max_active_power,
            )
            
            # Add to container using nested dictionary structure
            ts_container[resource_name][string(year)] = ts
        end
    end
    
    return ts_container
end

function create_ThermalStandard_PSY_timeseries(gen_variability_ts::DataFrame, ThermalStandard_collection::Vector{ThermalStandard})
    # Initialize container for storing timeseries as nested dictionary
    ts_container = Dict{String, Dict{String, SingleTimeSeries}}()
    
    # Define fixed calendar year timestamps (2035)
    tstamps = collect(range(DateTime("2035-01-01T00:00:00"), DateTime("2035-12-31T23:00:00"), step = Dates.Hour(1)))
    
    # Get unique years from the DateTime column and sort them
    years = sort(unique(Dates.year.(gen_variability_ts.DateTime)))
    
    # Loop through each renewable dispatch generator
    #for renewable in Renew_ND_collection[1:1] # testing with first element
    for thermal_standard in ThermalStandard_collection
        # Get the resource name and base power
        resource_name = get_name(thermal_standard)
        base_power = 1 # profiles are already normalized to base power
        
        # Initialize inner dictionary for this resource
        ts_container[resource_name] = Dict{String, SingleTimeSeries}()
        
        # Loop through each year in sorted order
        for year in years
            # Filter data for this year and sort by DateTime
            year_data = sort(gen_variability_ts[Dates.year.(gen_variability_ts.DateTime) .== year, :], :DateTime)
            
            # Get the availability column for this generator
            if !hasproperty(year_data, Symbol(resource_name))
                @warn "No availability data found for generator $resource_name in year $year"
                continue
            end
            
            # Verify hourly resolution, accounting for leap day skips
            for i in 2:length(year_data.DateTime)
                time_diff = year_data.DateTime[i] - year_data.DateTime[i-1]
                is_leap_day_skip = (Dates.month(year_data.DateTime[i-1]) == 2 && Dates.day(year_data.DateTime[i-1]) == 28 &&
                                  Dates.month(year_data.DateTime[i]) == 3 && Dates.day(year_data.DateTime[i]) == 1)
                
                if time_diff != Dates.Hour(1) && !is_leap_day_skip
                    @warn "Non-hourly resolution detected in year $year for $resource_name at index $i"
                end
            end
            
            # Normalize the data by base_power
            normalized_data = year_data[!, Symbol(resource_name)] ./ base_power
            
            # Create the timeseries name with year suffix
            ts_name = "max_active_power_$year"
            
            # Create the SingleTimeSeries object using fixed 2035 timestamps
            ts = SingleTimeSeries(;
                name = ts_name,
                data = TimeArray(tstamps, normalized_data),
                scaling_factor_multiplier = get_max_active_power,
            )
            
            # Add to container using nested dictionary structure
            ts_container[resource_name][string(year)] = ts
        end
    end
    
    return ts_container
end

function create_Hydro_PSY_timeseries(gen_variability_ts::DataFrame, Hydro_collection::Vector{HydroDispatch})
    # Initialize container for storing timeseries as nested dictionary
    ts_container = Dict{String, Dict{String, SingleTimeSeries}}()
    
    # Define fixed calendar year timestamps (2035)
    tstamps = collect(range(DateTime("2035-01-01T00:00:00"), DateTime("2035-12-31T23:00:00"), step = Dates.Hour(1)))
    
    # Get unique years from the DateTime column and sort them
    years = sort(unique(Dates.year.(gen_variability_ts.DateTime)))
    
    # Loop through each hydro generator
    for hydro in Hydro_collection
        # Get the resource name and base power
        resource_name = get_name(hydro)
        base_power = 1 # profiles are already normalized to base power
        
        # Initialize inner dictionary for this resource
        ts_container[resource_name] = Dict{String, SingleTimeSeries}()
        
        # Loop through each year in sorted order
        for year in years
            # Filter data for this year and sort by DateTime
            year_data = sort(gen_variability_ts[Dates.year.(gen_variability_ts.DateTime) .== year, :], :DateTime)
            
            # Get the availability column for this generator
            if !hasproperty(year_data, Symbol(resource_name))
                @warn "No availability data found for generator $resource_name in year $year"
                continue
            end
            
            # Verify hourly resolution, accounting for leap day skips
            for i in 2:length(year_data.DateTime)
                time_diff = year_data.DateTime[i] - year_data.DateTime[i-1]
                is_leap_day_skip = (Dates.month(year_data.DateTime[i-1]) == 2 && Dates.day(year_data.DateTime[i-1]) == 28 &&
                                  Dates.month(year_data.DateTime[i]) == 3 && Dates.day(year_data.DateTime[i]) == 1)
                
                if time_diff != Dates.Hour(1) && !is_leap_day_skip
                    @warn "Non-hourly resolution detected in year $year for $resource_name at index $i"
                end
            end
            
            # Normalize the data by base_power
            normalized_data = year_data[!, Symbol(resource_name)] ./ base_power
            
            # Create the timeseries name with year suffix
            ts_name = "max_active_power_$year"
            
            # Create the SingleTimeSeries object using fixed 2035 timestamps
            ts = SingleTimeSeries(;
                name = ts_name,
                data = TimeArray(tstamps, normalized_data),
                scaling_factor_multiplier = get_max_active_power,
            )
            
            # Add to container using nested dictionary structure
            ts_container[resource_name][string(year)] = ts
        end
    end
    
    return ts_container
end

function create_Hydro_Budget_PSY_timeseries(hydro_budget_ts::DataFrame, Hydro_collection::Vector{HydroDispatch})
    # Initialize container for storing timeseries as nested dictionary
    ts_container = Dict{String, Dict{String, SingleTimeSeries}}()
    
    # Define fixed calendar year timestamps (2035)
    tstamps = collect(range(DateTime("2035-01-01T00:00:00"), DateTime("2035-12-31T23:00:00"), step = Dates.Hour(1)))
    
    # Get unique years from the DateTime column and sort them
    years = sort(unique(Dates.year.(hydro_budget_ts.DateTime)))
    
    # Loop through each hydro generator
    for hydro in Hydro_collection
        # Get the resource name
        resource_name = get_name(hydro)
        
        # Initialize inner dictionary for this resource
        ts_container[resource_name] = Dict{String, SingleTimeSeries}()
        
        # Loop through each year in sorted order
        for year in years
            # Filter data for this year and sort by DateTime
            year_data = sort(hydro_budget_ts[Dates.year.(hydro_budget_ts.DateTime) .== year, :], :DateTime)
            
            # Get the budget column for this generator
            if !hasproperty(year_data, Symbol(resource_name))
                @warn "No budget data found for generator $resource_name in year $year"
                continue
            end
            
            # Get the budget data (already in p.u.)
            budget_data = year_data[!, Symbol(resource_name)]
            
            # Create the timeseries name with year suffix
            ts_name = "hydro_budget_$year"
            
            # Create the SingleTimeSeries object using fixed 2035 timestamps
            ts = SingleTimeSeries(;
                name = ts_name,
                data = TimeArray(tstamps, budget_data),
                scaling_factor_multiplier = get_max_active_power,
            )
            
            # Add to container using nested dictionary structure
            ts_container[resource_name][string(year)] = ts
        end
    end
    
    return ts_container
end


function create_reg_reserve_PSY_timeseries(
    demand_ts_df::DataFrame,
    gen_variability_ts_df::DataFrame,
    reserves_df::DataFrame,
    oprsv_zones::Vector{Int64},
    zone_dict::OrderedDict{String,Int64},
    Renew_D_generators::Vector{RenewableDispatch}
)
    # Initialize containers for storing timeseries with years as keys
    ts_container_up = Dict{String, SingleTimeSeries}()
    ts_container_down = Dict{String, SingleTimeSeries}()
    
    # Define fixed calendar year timestamps (2035)
    tstamps = collect(range(DateTime("2035-01-01T00:00:00"), DateTime("2035-12-31T23:00:00"), step = Dates.Hour(1)))
    
    # Get unique years from the DateTime column and sort them
    years = sort(unique(Dates.year.(demand_ts_df.DateTime)))
    
    # Initialize DataFrame to store the components
    reg_reserve_df = DataFrame(
        DateTime = DateTime[],
        Year = Int[],
        Demand_Component = Float64[],
        VRE_Component = Float64[],
        Total_Requirement = Float64[]
    )
    
    # Get the percentage requirements
    demand_percent = reserves_df[!, "Reg_Req_Percent_Demand"]
    vre_percent = reserves_df[!, "Reg_Req_Percent_VRE"]
    
    # Create reverse mapping from zone numbers to zone names
    zone_number_to_name = Dict{Int,String}()
    for (name, number) in zone_dict
        zone_number_to_name[number] = name
    end
    
    # Define the specific utility areas we want to include
    utility_areas = ["PGE", "SCE", "SDGE"]
    
    # Loop through each year
    for year in years
        # Filter data for this year
        year_demand = sort(demand_ts_df[Dates.year.(demand_ts_df.DateTime) .== year, :], :DateTime)
        year_gen = sort(gen_variability_ts_df[Dates.year.(gen_variability_ts_df.DateTime) .== year, :], :DateTime)
        
        # Initialize arrays for this year's components
        demand_component = zeros(length(tstamps))
        vre_component = zeros(length(tstamps))
        
        # Calculate demand-based component
        for zone_number in oprsv_zones
            # Convert zone number to zone name
            zone_name = zone_number_to_name[zone_number]
            demand_col = Symbol(zone_name * "_Demand")
            
            if hasproperty(year_demand, demand_col)
                @info "Adding demand data from Area:$zone_name to demand component of reserve requirement"
                demand_component .+= year_demand[!, demand_col] .* demand_percent
            else
                @warn "No demand data found for zone $zone_name (zone number $zone_number)"
            end
        end

        demand_component = round.(demand_component, digits=2)
        
        # Calculate VRE-based component
        # loop through all ACTIVE renewables in the system
        for gen in Renew_D_generators
            gen_name = get_name(gen)
            bus_name = get_name(get_bus(gen))
            
            # Only include VRE in the specified utility areas
            if bus_name in utility_areas && hasproperty(year_gen, Symbol(gen_name))
                @info "Adding VRE data from generator:$gen_name in Area:$bus_name to VRE component of reserve requirement"
                vre_component .+= year_gen[!, Symbol(gen_name)] .* get_base_power(gen) .* vre_percent
            end
        end

        vre_component = round.(vre_component, digits=2)
        
        # Calculate total requirement
        total_requirement = demand_component .+ vre_component
        total_requirement = round.(total_requirement, digits=2)

        # Add to components DataFrame
        for i in 1:length(tstamps)
            push!(reg_reserve_df, (
                tstamps[i],
                year,
                demand_component[i],
                vre_component[i],
                total_requirement[i]
            ))
        end
        
        # Create the SingleTimeSeries objects for up and down reserves (each with half the requirement)
        ts_up = SingleTimeSeries(
            name = "requirement_up_$year",
            data = TimeArray(tstamps, total_requirement ./ 2), # dividing by 2 to get half the requirement
            scaling_factor_multiplier = get_requirement, # defining in natural units
        )
        
        ts_down = SingleTimeSeries(
            name = "requirement_down_$year",
            data = TimeArray(tstamps, total_requirement ./ 2), # dividing by 2 to get half the requirement
            scaling_factor_multiplier = get_requirement, # defining in natural units
        )
        
        # Add to respective containers using year as key
        ts_container_up[string(year)] = ts_up
        ts_container_down[string(year)] = ts_down
    end
    
    return ts_container_up, ts_container_down, reg_reserve_df
end

function create_fuel_price_PSY_timeseries(fuel_ts::DataFrame)
    # Initialize container for storing timeseries by fuel type
    ts_container = Dict{String, SingleTimeSeries}()
    
    # Define fixed calendar year timestamps (2035)
    tstamps = collect(range(DateTime("2035-01-01T00:00:00"), DateTime("2035-12-31T23:00:00"), step = Dates.Hour(1)))
    
    # Get all fuel type columns from the DataFrame
    fuel_columns = names(fuel_ts, Not(:DateTime))
    
    # Filter data for first year only (1998)
    first_year_data = fuel_ts[Dates.year.(fuel_ts.DateTime) .== 1998, :]
    
    # Loop through each fuel type
    for fuel_col in fuel_columns
        # Get the fuel price data for first year only
        fuel_price_data = first_year_data[!, fuel_col]
        
        # Create the timeseries name
        ts_name = "fuel_price"
        
        # Create the SingleTimeSeries object using fixed 2035 timestamps
        ts = SingleTimeSeries(;
            name = ts_name,
            data = TimeArray(tstamps, fuel_price_data),
            #scaling_factor_multiplier = 1.0,
        )
        
        # Add to container using fuel type as key
        ts_container[string(fuel_col)] = ts
    end

    #ts_container["SW_Coal_Fuel"]
    
    return ts_container
end

function create_generic_requirement_reserveUp_timeseries(sys::System, WY::Int64)
    # first we need to remove all forecasts (i.e. DeterministicSingleTimeSeries) from the system
    #remove_time_series!(sys, DeterministicSingleTimeSeries)
    
    # Get the reserve up service
    reserve_up = get_component(VariableReserve{ReserveUp}, sys, "CAISO_reg_up")
    
    # Check if "requirement" timeseries already exists and remove it if it does
    if "requirement" ∈ get_name.(get_time_series_keys(reserve_up))
        remove_time_series!(sys, SingleTimeSeries, reserve_up, "requirement")
    else
        # do nothing
    end
    
    # Get the timeseries array for the specific weather year
    ts_array = get_time_series_array(SingleTimeSeries, reserve_up, "requirement_up_$WY"; ignore_scaling_factors = true)
    
    # Create a new timeseries with the same data but named "requirement"
    tstamp = timestamp(ts_array)
    vals = values(ts_array)
    new_ts = SingleTimeSeries(
        name = "requirement",
        data = TimeArray(tstamp, vals),
        scaling_factor_multiplier = get_requirement
    )
    
    # Add the new timeseries to the reserve up service
    add_time_series!(sys, reserve_up, new_ts)    
end

function create_generic_requirement_reserveDown_timeseries(sys::System, WY::Int64)
    # first we need to remove all forecasts (i.e. DeterministicSingleTimeSeries) from the system
    #remove_time_series!(sys, DeterministicSingleTimeSeries)
    
    # Get the reserve down service
    reserve_down = get_component(VariableReserve{ReserveDown}, sys, "CAISO_reg_down")
    
    # Check if "requirement" timeseries already exists and remove it if it does
    if "requirement" ∈ get_name.(get_time_series_keys(reserve_down))
        remove_time_series!(sys, SingleTimeSeries, reserve_down, "requirement")
    else
        # do nothing
    end
    
    # Get the timeseries array for the specific weather year
    ts_array = get_time_series_array(SingleTimeSeries, reserve_down, "requirement_down_$WY"; ignore_scaling_factors = true)
    
    # Create a new timeseries with the same data but named "requirement"
    tstamp = timestamp(ts_array)
    vals = values(ts_array)
    new_ts = SingleTimeSeries(
        name = "requirement",
        data = TimeArray(tstamp, vals),
        scaling_factor_multiplier = get_requirement
    )
    
    # Add the new timeseries to the reserve down service
    add_time_series!(sys, reserve_down, new_ts)    
end

function create_generic_daily_hydrobudget_timeseries(sys::System, WY::Int64, hydro_collection::Vector{HydroDispatch})
    # first we need to remove all forecasts (i.e. DeterministicSingleTimeSeries) from the system
    #remove_time_series!(sys, DeterministicSingleTimeSeries)

    # Loop through each hydro unit
    for hydro in hydro_collection
        # Check if "hydro_budget" timeseries already exists and remove it if it does
        if "hydro_budget" ∈ get_name.(get_time_series_keys(hydro))
            remove_time_series!(sys, SingleTimeSeries, hydro, "hydro_budget")
        end
        
        # Get the timeseries array for the specific weather year
        ts_array = get_time_series_array(SingleTimeSeries, hydro, "hydro_budget_$WY"; ignore_scaling_factors = true)
        
        # Create a new timeseries with the same data but named "hydro_budget"
        tstamp = timestamp(ts_array)
        vals = values(ts_array)
        new_ts = SingleTimeSeries(
            name = "hydro_budget",
            data = TimeArray(tstamp, vals),
            scaling_factor_multiplier = get_max_active_power
        )
        
        # Add the new timeseries to the hydro unit
        add_time_series!(sys, hydro, new_ts)
    end
end


function create_PHS_PSY_timeseries(PHS_collection::Vector{HydroPumpedStorage})
    # Initialize container for storing timeseries as nested dictionary
    ts_container = Dict{String, Dict{String, SingleTimeSeries}}()
    
    # Define fixed calendar year timestamps (2035)
    tstamps = collect(range(DateTime("2035-01-01T00:00:00"), DateTime("2035-12-31T23:00:00"), step = Dates.Hour(1)))
    
    # Create zero-valued array for all timestamps
    zero_data = zeros(length(tstamps))
    
    # Loop through each PHS unit
    for phs in PHS_collection
        # Get the resource name
        resource_name = get_name(phs)
        
        # Initialize inner dictionary for this resource
        ts_container[resource_name] = Dict{String, SingleTimeSeries}()
        
        # Create inflow timeseries
        inflow_ts = SingleTimeSeries(
            name = "inflow",
            data = TimeArray(tstamps, zero_data),
            scaling_factor_multiplier = get_inflow,
        )
        
        # Create outflow timeseries
        outflow_ts = SingleTimeSeries(
            name = "outflow",
            data = TimeArray(tstamps, zero_data),
            scaling_factor_multiplier = get_outflow,
        )
        
        # Add to container using nested dictionary structure
        ts_container[resource_name]["inflow"] = inflow_ts
        ts_container[resource_name]["outflow"] = outflow_ts
    end
    
    return ts_container
end

function update_TS_fuel_price!(sys::System, thermal_standards::Vector{ThermalStandard})
    # reassign fuel_price timeseries to ThermalStandard objects
    for g in thermal_standards
        if "fuel_price" ∈ get_name.(get_time_series_keys(g))
            # fuel_ts = get_time_series(SingleTimeSeries, g, "fuel_price")
            fuel_array = get_time_series_array(SingleTimeSeries, g, "fuel_price"; ignore_scaling_factors = true)
            tstamp = timestamp(fuel_array)
            vals = values(fuel_array)
            new_ts = SingleTimeSeries("fuel_price", TimeArray(tstamp, vals))
            remove_time_series!(sys, SingleTimeSeries, g, "fuel_price")
            set_fuel_cost!(sys, g, new_ts)
        end
    end
end

function update_TS_PHS_flows!(sys::System, PHS_objects::Vector{HydroPumpedStorage})
    # reassign inflow and outflow timeseries to PHS objects
    for phs in PHS_objects
        # Handle inflow timeseries
        if "inflow" ∈ get_name.(get_time_series_keys(phs))
            inflow_array = get_time_series_array(SingleTimeSeries, phs, "inflow"; ignore_scaling_factors = true)
            tstamp = timestamp(inflow_array)
            vals = values(inflow_array)
            new_inflow_ts = SingleTimeSeries("inflow", TimeArray(tstamp, vals))
            remove_time_series!(sys, SingleTimeSeries, phs, "inflow")
            set_inflow!(sys, phs, new_inflow_ts)
        end

        # Handle outflow timeseries
        if "outflow" ∈ get_name.(get_time_series_keys(phs))
            outflow_array = get_time_series_array(SingleTimeSeries, phs, "outflow"; ignore_scaling_factors = true)
            tstamp = timestamp(outflow_array)
            vals = values(outflow_array)
            new_outflow_ts = SingleTimeSeries("outflow", TimeArray(tstamp, vals))
            remove_time_series!(sys, SingleTimeSeries, phs, "outflow")
            set_outflow!(sys, phs, new_outflow_ts)
        end
    end
end



