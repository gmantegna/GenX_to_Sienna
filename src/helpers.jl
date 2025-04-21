function initialize_paths_and_inputs()
    # Define root directory and related paths
    current_dir = pwd()

    # Check if "Sonoma" exists in the path (want to have root directory get you to Sonoma)
    if occursin("Sonoma", current_dir)
        # Find the position where "Sonoma" starts
        start_pos = first(findfirst("Sonoma", current_dir))
        # Trim off everything from "Sonoma" onward
        root_dir = current_dir[1:start_pos-1]
    else
        root_dir = current_dir
    end

    # Define key directories
    data_dir = joinpath(root_dir, "Sonoma", "GENX_Output") # GenX input files
    output_dir_base = joinpath(root_dir, "Sonoma", "Sienna_Outputs") # Base directory for outputs
    
    # Ensure output directory exists
    if !ispath(output_dir_base)
        mkpath(output_dir_base)
    end

    # Return all paths as a dictionary for easy access
    return Dict(
        :root_dir => root_dir,
        :data_dir => data_dir,
        :output_dir_base => output_dir_base,
        # :genx_capacity_results => joinpath(data_dir, "results")  # This is a constant from the original code
    )
end

function process_demand_data(demand_data_path::String, zones_dict::OrderedDict{String,Int64})
    # Read the demand data
    demand_df = CSV.read(demand_data_path, DataFrame)

    # Create a mapping from zone numbers to zone names
    zone_number_to_name = Dict{Int,String}()
    for (name, number) in zones_dict
        zone_number_to_name[number] = name
    end
    
    # only keep the demand columns (i.e. prefix with "Demand_MW_z")
    # retrieve column for load zones and convert to string
    demand_columns = names(demand_df)[startswith.(string.(names(demand_df)), "Demand_MW_z")]
    
    # Create new column names mapping
    new_column_names = Dict{String,String}()
    for col in demand_columns
        # Extract zone number from column name (e.g., "Demand_MW_z1" -> 1)
        zone_num = parse(Int, match(r"Demand_MW_z(\d+)", col).captures[1])
        # Get the actual zone name from our mapping
        zone_name = zone_number_to_name[zone_num]
        # Create new column name
        new_column_names[col] = zone_name * "_Demand"
    end
    
    # select only the demand columns and rename them
    select!(demand_df, demand_columns)
    DataFrames.rename!(demand_df, new_column_names)
    
    # Create datetime vector for 1998-2020 (23 years)
    dates = Vector{DateTime}()
    demand_zones = Dict{String, Vector{Float64}}()
    
    # Initialize demand vectors for each zone using the new column names
    for col_name in names(demand_df)
        demand_zones[String(col_name)] = Float64[]
    end
    
    # Process each year
    for year in 1998:2020
        # Determine number of days in year
        days_in_year = Dates.isleapyear(year) ? 366 : 365
        
        for day in 1:days_in_year
            # Skip leap days
            if Dates.isleapyear(year) && day == 60  # February 29
                continue
            end
            
            # Add 24 hours for this day
            for hour in 0:23
                current_date = DateTime(year, 1, 1) + Dates.Day(day-1) + Dates.Hour(hour)
                push!(dates, current_date)
                
                # Calculate index in the original data
                # Need to account for skipped leap days in previous years
                leap_days_before = sum(Dates.isleapyear.(1998:year-1))
                if Dates.isleapyear(year) && day > 59  # After Feb 28 in leap year
                    leap_days_before += 1
                end
                
                data_index = ((year - 1998) * 365 + day - 1) * 24 + hour + 1 - leap_days_before
                
                # For 12/31/2020, use 12/30/2020 data
                if year == 2020 && day == 365
                    data_index = data_index - 24  # Use previous day's data
                end
                
                # Add demand for each zone using new column names
                for col_name in names(demand_df)
                    push!(demand_zones[String(col_name)], demand_df[data_index, col_name])
                end
            end
        end
    end
    
    # Create final DataFrame
    result_df = DataFrame(:DateTime => dates) # initialize the DataFrame with the dates
    for (zone_name, values) in demand_zones # loop through the demand zones and add the values to the DataFrame
        result_df[!, zone_name] = values
    end
    
    return result_df
end

function process_fuel_data(fuels_df::DataFrame)
    
    # Drop the first row (CO2 emissions) and the first column (index)
    fuels_df = fuels_df[2:end, 2:end]
    
    # Get all column names except the first (which is likely an index)
    fuel_columns = names(fuels_df)
    
    # Create datetime vector for 1998-2020 (23 years)
    dates = Vector{DateTime}()
    fuel_prices = Dict{String, Vector{Float64}}()
    
    # Initialize fuel price vectors using the actual column names
    for col_name in fuel_columns
        fuel_prices[String(col_name)] = Float64[]
    end
    
    # Process each year
    for year in 1998:2020
        # Determine number of days in year
        days_in_year = Dates.isleapyear(year) ? 366 : 365
        
        for day in 1:days_in_year
            # Skip leap days
            if Dates.isleapyear(year) && day == 60  # February 29
                continue
            end
            
            # Add 24 hours for this day
            for hour in 0:23
                current_date = DateTime(year, 1, 1) + Dates.Day(day-1) + Dates.Hour(hour)
                push!(dates, current_date)
                
                # Calculate index in the original data
                # Need to account for skipped leap days in previous years
                leap_days_before = sum(Dates.isleapyear.(1998:year-1))
                if Dates.isleapyear(year) && day > 59  # After Feb 28 in leap year
                    leap_days_before += 1
                end
                
                data_index = ((year - 1998) * 365 + day - 1) * 24 + hour + 1 - leap_days_before
                
                # For 12/31/2020, use 12/30/2020 data
                if year == 2020 && day == 365
                    data_index = data_index - 24  # Use previous day's data
                end
                
                # Add fuel prices using actual column names
                for col_name in fuel_columns
                    push!(fuel_prices[String(col_name)], fuels_df[data_index, col_name])
                end
            end
        end
    end
    
    # Create final DataFrame
    result_df = DataFrame(:DateTime => dates)
    for (fuel_name, values) in fuel_prices
        result_df[!, fuel_name] = values
    end
    
    return result_df
end

function process_generator_variability_data(gen_var_data_path::String)
    # Read the generator variability data
    gen_var_df = CSV.read(gen_var_data_path, DataFrame)
    
    # Drop the first column (time index)
    select!(gen_var_df, Not(1))
    
    # Get all generator columns
    generator_columns = names(gen_var_df)
    
    # Create datetime vector for 1998-2020 (23 years)
    dates = Vector{DateTime}()
    generator_profiles = Dict{String, Vector{Float64}}()
    
    # Initialize generator profile vectors using the actual column names
    for col_name in generator_columns
        generator_profiles[String(col_name)] = Float64[]
    end
    
    # Process each year
    for year in 1998:2020
        # Determine number of days in year
        days_in_year = Dates.isleapyear(year) ? 366 : 365
        
        for day in 1:days_in_year
            # Skip leap days
            if Dates.isleapyear(year) && day == 60  # February 29
                continue
            end
            
            # Add 24 hours for this day
            for hour in 0:23
                current_date = DateTime(year, 1, 1) + Dates.Day(day-1) + Dates.Hour(hour)
                push!(dates, current_date)
                
                # Calculate index in the original data
                # Need to account for skipped leap days in previous years
                leap_days_before = sum(Dates.isleapyear.(1998:year-1))
                if Dates.isleapyear(year) && day > 59  # After Feb 28 in leap year
                    leap_days_before += 1
                end
                
                data_index = ((year - 1998) * 365 + day - 1) * 24 + hour + 1 - leap_days_before
                
                # For 12/31/2020, use 12/30/2020 data
                if year == 2020 && day == 365
                    data_index = data_index - 24  # Use previous day's data
                end
                
                # Add generator profiles using actual column names
                for col_name in generator_columns
                    push!(generator_profiles[String(col_name)], gen_var_df[data_index, col_name])
                end
            end
        end
    end
    
    # Create final DataFrame
    result_df = DataFrame(:DateTime => dates)
    for (generator_name, values) in generator_profiles
        result_df[!, generator_name] = values
    end
    
    return result_df
end

