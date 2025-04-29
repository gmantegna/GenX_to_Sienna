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


