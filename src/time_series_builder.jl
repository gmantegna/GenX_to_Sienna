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