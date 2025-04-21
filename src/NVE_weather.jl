function define_FTM_solar_time_series(solar_dir::String, sys::System)
    #################
    # Step 0: Check units base
    #################
    get_units_base(sys)
    set_units_base_system!(sys, "NATURAL_UNITS")

    #################
    # Step 1: Read In Solar Data
    #################
    # Create a dictionary to store the data, organizing it by year
    ftm_solar_ts_dict = Dict{String, DataFrame}()

    # process all .csv files in the specified directory
    for file in readdir(solar_dir)
        if endswith(file, ".csv")
            # Extract the year from the filename (e.g., "Y1998.csv" -> "Y1998")
            wy = splitext(file)[1]

            # Read the file into a DataFrame
            file_path = joinpath(solar_dir, file)
            df = CSV.read(file_path, DataFrame)

            # Filter out leap day (February 29)
            df = filter(row -> !(row.Month == 2 && row.Day == 29), df)

            # Store the filtered DataFrame in the dictionary
            ftm_solar_ts_dict[wy] = df
        end
    end

    ##########################
    # Step 2: Pull Up RenewableDispatch Solar Resources
    ##########################
    # Filter solar RenewableDispatch resources from system (this returns an iterable object)
    #solar_resources = get_components(x -> get_prime_mover_type(x) == PrimeMovers.PVe, RenewableGen, sys)
    FTM_solar_resources = get_components(x -> x isa RenewableDispatch && get_prime_mover_type(x) == PrimeMovers.PVe, RenewableGen, sys)

    # Define an 8760 vector based on a date interval object
    tstamps = collect(range(DateTime("2030-01-01T00:00:00"), DateTime("2030-12-31T23:00:00"), step = Dates.Hour(1)))

    # MANUAL INPUT: Solar mapping for missing profiles (active_object: mapped_object)
    solar_device_mapping = Dict(
        "Dry Lake East Solar" => "Dry Lake Solar",
        "_NVE IRP Expansion_SPPC_PV (Standalone)" => "NV Generic Solar Hybrid",
        "_NVE IRP Expansion_NVP_PV (Hybrid)" => "NV Generic Solar Hybrid",
        "_PCM Generic Expansion_NVP_PV" => "NV Generic Solar",
        "_PCM Generic Expansion_SPPC_PV" => "NV Generic Solar",
        "Libra Solar" => "NV Generic Solar",
        "Sierra Solar" => "NV Generic Solar Hybrid"        
    )

    ###########################
    # Step 3: Retrieve the list of generator names in the FTM CSV folder
    ###########################
    active_df = first(values(ftm_solar_ts_dict))  # grab the first DataFrame in dictionary to use as reference
    ftm_solar_ts_dict_s = String[]  # Initialize an array to store names of active solar facilities defined in the csvs

    # Define ftm_solar_ts_dict_s by looping through all column names in the DataFrame
    for name in names(active_df)
        if name ∈ ["Month", "Day", "Period"]
            continue  # Skip these column names
        else
            push!(ftm_solar_ts_dict_s, name)  # Append valid names
        end
    end

    ##############################
    # Step 4: Record Annual Max Hourly CF value across all weather years
    ##############################
    # Initialize a DataFrame to store results
    df_max_solar = DataFrame(
        FTM_solar_resource = String[],
        max_pu_value = Float64[],
        weather_yr = String[])

    for device in FTM_solar_resources #loop through FTM solar resources
        # Retrieve device name
        device_name = get_name(device)

        # initial placeholder values
        max_value = -Inf
        max_year = ""

        # Reset loop variables
        selected_name = nothing

        # Define the profile to assign to the device
        if device_name in ftm_solar_ts_dict_s
            selected_name = device_name
        elseif haskey(solar_device_mapping, device_name)
            substituted_name = solar_device_mapping[device_name]
            if substituted_name in ftm_solar_ts_dict_s
                selected_name = substituted_name
            end
        end

        # If still no match found after applying the manual mappings, assign default time series and log a warning
        if selected_name === nothing
            selected_name = "NV Generic Solar"
            @warn "No time series found for device $device_name. Using default time series: $selected_name."
        end

        # Collect all time series for the active device by looping through all weather years
        for (wy, data) in ftm_solar_ts_dict
            # println("processing weather year: ", wy)

            #testing / debugging
            # wy = "Y1998"
            # data = solar_ts_dict[wy]

            if selected_name in names(data)  # Check if the selected name exists in the DataFrame

                # pull max value of current unit for current wy
                local_max = maximum(data[!, Symbol(selected_name)])

                if local_max > max_value
                    # println("new max value; updating df")
                    max_value = local_max
                    max_year = wy
                end
            else
                # log a warning if something is wrong
                @warn "Selected name $selected_name not found in data for year $year."
            end
        end #looping through weather years
        # Append the results to the DataFrame
        push!(df_max_solar, (device_name, max_value, max_year))
    end # looping through solar resources

    #df_max_solar

    ##############################
    # Step 5: create and store PowerSystems.jl timeseries
    ##############################
    # initialize container
    FTM_ts_container = Vector{Tuple{String, SingleTimeSeries}}()  # Container for all time series of all FTM solar devices

    for device in FTM_solar_resources #loop through solar resources

        #testing /debugging
        #device = collect(FTM_solar_resources)[1]

        # Retrieve device name
        device_name = get_name(device)

        # Reset loop variables
        selected_name = nothing

        # Define the profile to assign to the device
        if device_name in ftm_solar_ts_dict_s  #is the active device in the csvs
            selected_name = device_name
        elseif haskey(solar_device_mapping, device_name) #otherwise use the manually defined mapping
            substituted_name = solar_device_mapping[device_name]
            if substituted_name in ftm_solar_ts_dict_s
                selected_name = substituted_name
            end
        end

        # If still no match found, assign default time series and log a warning
        if selected_name === nothing
            selected_name = "NV Generic Solar"
            # @warn "No time series found for device $device_name. Using default time series: $selected_name."
        end

        # define max annual hourly cf across all weather years (mininum of 1)
        max_cf_all_wys = max(1,df_max_solar[df_max_solar.FTM_solar_resource .== device_name, :max_pu_value][1])
        # retrieve max_active_power (default weather year)
        max_ap = get_max_active_power(device)
        # calculate revised max_active_power
        #max_ap_rev = max_ap * max_cf_all_wys
        # reassign rating (i.e., max_active_power) parameter
        #set_rating!(device, max_ap_rev)

        # Collect all time series for the active device by looping through all weather years
        for (wy, data) in ftm_solar_ts_dict
            # println("processing weather year: ", wy)

            #testing / debugging
            #wy = "Y1998"
            #data = ftm_solar_ts_dict[wy]

            # Generate a unique time series name
            ts_name = "max_active_power_$(wy)"  # Unique name for each time series

            # normalizing data to reflect revised max_active_power (i.e. rating)
            normalized_data = data[!, selected_name] # raw data is provided in p.u. units

            if selected_name in names(data)  # Check if the selected name exists in the DataFrame
                # Create the time series by extracting data from the DataFrame
                ts = SingleTimeSeries(;
                    name = ts_name,
                    data = TimeArray(tstamps, normalized_data), #create a TimeArray object w/ p.u. values
                    scaling_factor_multiplier = get_max_active_power, #scaling factor multiplier (we need this for SiennaPRASInterface runs)
                )
                # Store the time series along with its device name
                push!(FTM_ts_container, (device_name, ts))
            else
                # Log a warning if something is wrong
                @warn "Selected name $selected_name not found in data for year $year."
            end
        end # looping through weather years
    end # looping through solar resources
    return FTM_ts_container, df_max_solar
end

function add_FTM_solar_time_series_to_system!(sys::System, FTM_s_ts_container::Vector{Tuple{String, SingleTimeSeries}})
    # Loop through each solar RenewableDispatch object and add time series
    for (device_name, time_series) in FTM_s_ts_container
        # Retrieve the active device by its name
        active_device = get_component(RenewableDispatch, sys, device_name)

        if active_device !== nothing
            # Add the time series to the system
            add_time_series!(sys, active_device, time_series)
            # println("Added time series: ", time_series.name, " to device: ", device_name)
        else
            @warn "Device $device_name not found in the system. Time series $(time_series.name) not added."
        end
    end
end

##################################################################################

function define_BTM_solar_time_series(solar_dir::String, sys::System)
    #################
    # Step 0: Check units base
    #################
    get_units_base(sys)
    set_units_base_system!(sys, "NATURAL_UNITS")

    #################
    # Step 1: Read In BTM solar data - which is in 30 min increments and convert to hourly
    #################
    # Create a dictionary to store the data, organizing it by year
    btm_solar_ts_dict = Dict{String, DataFrame}()

    # process all .csv files in the specified directory
    for file in readdir(solar_dir)
        if endswith(file, ".csv")
            # Extract the year from the filename (e.g., "Y1998.csv" -> "Y1998")
            wy = splitext(file)[1]

            # Read the file into a DataFrame
            file_path = joinpath(solar_dir, file)
            df = CSV.read(file_path, DataFrame)

            # Filter out leap day (February 29)
            df = filter(row -> !(row.Month == 2 && row.Day == 29), df)

            # Add Hour Ending (HE) column: maps Period 1–2 → HE 1, 3–4 → HE 2, ..., 47–48 → HE 24
            df.HE = ceil.(Int, df.Period ./ 2)

            # Group by Month, Day, and HE and compute the hourly average
            df = combine(groupby(df, [:Month, :Day, :HE]), :DPV_NV => mean => :DPV_NV)

            # rename the HE column to Period
            DataFrames.rename!(df, :HE => :Period)

            # Store the filtered DataFrame in the dictionary
            btm_solar_ts_dict[wy] = df
        end
    end

    ##########################
    # Step 2: Pull Up RenewableNONDispatch Solar Resources
    ##########################
    # Filter for RenewableNon Dispatch for solar (this returns an iterable object)
    btm_solar_resources = get_components(x -> x isa RenewableNonDispatch && get_prime_mover_type(x) == PrimeMovers.PVe, RenewableGen, sys)

    # Define an 8760 vector based on a date interval object
    tstamps = collect(range(DateTime("2030-01-01T00:00:00"), DateTime("2030-12-31T23:00:00"), step = Dates.Hour(1)))

    ###########################
    # Step 3: Retrieve the list of generator names in the btm CSV folder
    ###########################
    active_df = first(values(btm_solar_ts_dict))  # grab the first DataFrame in dictionary to use as reference
    btm_solar_ts_dict_s = String[]  # Initialize an array to store names of active solar facilities defined in the csvs

    # Define btm_solar_ts_dict_s by looping through all column names in the DataFrame
    for name in names(active_df)
        if name ∈ ["Month", "Day", "Period"]
            continue  # Skip these column names
        else
            push!(btm_solar_ts_dict_s, name)  # Append valid names
        end
    end

    ##############################
    # Step 4: Record Annual Max Hourly CF value across all weather years
    ##############################
    # Initialize a DataFrame to store results
    df_max_BTM_solar = DataFrame(
        btm_solar_resource = String[],
        max_pu_value = Float64[],
        weather_yr = String[])

    for device in btm_solar_resources #loop through btm solar resources
        # Retrieve device name
        device_name = get_name(device)

        # initial placeholder values
        max_value = -Inf
        max_year = ""

        # Reset loop variables
        selected_name = nothing

        # Define the profile to assign to the device
        if device_name in btm_solar_ts_dict_s
            selected_name = device_name
        else
            selected_name = "DPV_NV" #default name for btm solar resources
        end

        # Collect all time series for the active device by looping through all weather years
        for (wy, data) in btm_solar_ts_dict
            # println("processing weather year: ", wy)

            #testing / debugging
            # wy = "Y1998"
            # data = solar_ts_dict[wy]

            if selected_name in names(data)  # Check if the selected name exists in the DataFrame

                # pull max value of current unit for current wy
                local_max = maximum(data[!, Symbol(selected_name)])

                if local_max > max_value
                    # println("new max value; updating df")
                    max_value = local_max
                    max_year = wy
                end
            else
                # log a warning if something is wrong
                @warn "Selected name $selected_name not found in data for year $year."
            end
        end #looping through weather years
        # Append the results to the DataFrame
        push!(df_max_BTM_solar, (device_name, max_value, max_year))
    end # looping through solar resources

    #df_max_BTM_solar

    ##############################
    # Step 5: create and store PowerSystems.jl timeseries
    ##############################
    # initialize container
    btm_ts_container = Vector{Tuple{String, SingleTimeSeries}}()  # Container for all time series of all solar devices

    for device in btm_solar_resources #loop through solar resources

        # Retrieve device name
        device_name = get_name(device)

        # Reset loop variables
        selected_name = nothing

        # Define the profile to assign to the device
        if device_name in btm_solar_ts_dict_s  #is the active device in the csvs
            selected_name = device_name
        else
            selected_name = "DPV_NV" # default name for BTM PV
        end

        # define max annual hourly cf across all weather years (mininum of 1)
        max_cf_all_wys = max(1,df_max_BTM_solar[df_max_BTM_solar.btm_solar_resource .== device_name, :max_pu_value][1])
        # retrieve max_active_power (default weather year)
        max_ap = get_max_active_power(device)
        # calculate revised max_active_power
        #max_ap_rev = max_ap * max_cf_all_wys
        # reassign rating (i.e., max_active_power) parameter
        #set_rating!(device, max_ap_rev)

        # Collect all time series for the active device by looping through all weather years
        for (wy, data) in btm_solar_ts_dict
            # println("processing weather year: ", wy)

            #testing / debugging
            #wy = "Y1998"
            #data = btm_solar_ts_dict[wy]

            # Generate a unique time series name
            ts_name = "max_active_power_$(wy)"  # Unique name for each time series

            # normalizing data to reflect revised max_active_power (i.e. rating)
            normalized_data = data[!, selected_name] # raw data is provided in p.u. units

            if selected_name in names(data)  # Check if the selected name exists in the DataFrame
                # Create the time series by extracting data from the DataFrame
                ts = SingleTimeSeries(;
                    name = ts_name,
                    data = TimeArray(tstamps, normalized_data),#create a TimeArray object w/ p.u. values
                    scaling_factor_multiplier = get_max_active_power,#apply scaling_factor_multiplier (we need this for SiennaPRASInterface runs)
                )
                
                # Store the time series along with its device name
                push!(btm_ts_container, (device_name, ts))
            else
                # Log a warning if something is wrong
                @warn "Selected name $selected_name not found in data for year $year."
            end
        end # looping through weather years
    end # looping through solar resources
    return btm_ts_container, df_max_BTM_solar
end

function add_BTM_solar_time_series_to_system!(sys::System, btm_s_ts_container::Vector{Tuple{String, SingleTimeSeries}})
    # Loop through each solar RenewableDispatch object and add time series
    for (device_name, time_series) in btm_s_ts_container
        # Retrieve the active device by its name
        active_device = get_component(RenewableNonDispatch, sys, device_name)

        if active_device !== nothing
            # Add the time series to the system
            add_time_series!(sys, active_device, time_series)
            # println("Added time series: ", time_series.name, " to device: ", device_name)
        else
            @warn "Device $device_name not found in the system. Time series $(time_series.name) not added."
        end
    end
end

#####################################################################################

function define_wind_ts_container(wind_dir::String, sys::System)
    #################
    # Step 1: Wind Weather Year Mapping
    #################
    # Define the wind weather year mapping table (one-to-many assignment)
    wind_weather_mapping = Dict(
        "Y1998" => "2014",
        "Y1999" => "2007",
        "Y2000" => "2008",
        "Y2001" => "2009",
        "Y2002" => "2010",
        "Y2003" => "2011",
        "Y2004" => "2012",
        "Y2005" => "2013",
        "Y2006" => "2014",
        "Y2007" => "2007",
        "Y2008" => "2008",
        "Y2009" => "2009",
        "Y2010" => "2010",
        "Y2011" => "2011",
        "Y2012" => "2012",
        "Y2013" => "2013",
        "Y2014" => "2014",
        "Y2015" => "2007",
        "Y2016" => "2008",
        "Y2017" => "2009",
        "Y2018" => "2010",
        "Y2019" => "2011",
        "Y2020" => "2012",
        "Y2021" => "2013",
        "Y2022" => "2014"
    )

    #################
    # Step 2: Load Wind Data
    #################
    # Create a dictionary to store the original wind weather year data
    base_wind_ts_dict = Dict{String, DataFrame}()

    # Process all .csv files in the specified directory
    for file in readdir(wind_dir)
        if endswith(file, ".csv")
            # Extract the year from the filename (e.g., "2007.csv" -> "2007")
           wy = splitext(file)[1]

            # Read the file into a DataFrame
            file_path = joinpath(wind_dir, file)
            df = CSV.read(file_path, DataFrame)

            # Filter out leap day (February 29)
            df = filter(row -> !(row.Month == 2 && row.Day == 29), df)

            # Store the filtered DataFrame in the base dictionary
            base_wind_ts_dict[wy] = df
        end
    end

    # Create the final dictionary with 25 entries
    wind_ts_dict = Dict{String, DataFrame}()

    for (sienna_wy, wind_yr) in wind_weather_mapping
        # Use the mapping from wind_weather_mapping dictionary
        wind_key = "Y$wind_yr"
        if haskey(base_wind_ts_dict, wind_key)
            wind_ts_dict[sienna_wy] = base_wind_ts_dict[wind_key]
        else
            @warn "Wind weather year $wind_key not found. Sienna year $sienna_wy will not have data."
        end
    end

    #################
    # Step 3: Prepare Generator Names and Resources
    #################
    # Filter RenewableDispatch for wind
    wind_resources = get_components(x -> get_prime_mover_type(x) == PrimeMovers.WT, RenewableGen, sys)

    # Define an 8760 vector based on a date interval object
    tstamps = collect(range(DateTime("2030-01-01T00:00:00"), DateTime("2030-12-31T23:00:00"), step = Dates.Hour(1)))

    # MANUAL UPDATE: Wind mapping for missing profiles
    wind_device_mapping = Dict(
        "Idaho Wind" => "Wind_Idaho", #Wind_Idaho is the name of the column in the weather yr
        "_PCM Generic Expansion_SPPC_Wind (ID)" => "Wind_Idaho",
    )

    # Define the list of generator names in the CSV folder
    active_df = first(values(wind_ts_dict))  # Get the first DataFrame
    wind_ts_dict_s = String[]  # Initialize an array to store column names

    for name in names(active_df)
        if name ∈ ["Month", "Day", "Period"]
            continue  # Skip these column names
        else
            push!(wind_ts_dict_s, name)  # Append valid names
        end
    end

    #################
    # Step 4: Create Time Series Container
    #################
    ts_container = Vector{Tuple{String, SingleTimeSeries}}()  # Initialize container for wind time series

    for device in wind_resources
        # Retrieve device name
        device_name = get_name(device)
        # retrieve max_active_power
        max_power = get_max_active_power(device)

        # Reset loop variables
        selected_name = nothing

        # Define the profile to assign to the device
        if device_name in wind_ts_dict_s
            selected_name = device_name
        elseif haskey(wind_device_mapping, device_name)
            substituted_name = wind_device_mapping[device_name]
            if substituted_name in wind_ts_dict_s
                selected_name = substituted_name
            end
        end

        # If no match found, assign default time series and log the assignment
        if selected_name === nothing
            selected_name = "Wind_NV-AZ"
            @warn "No time series found for device $device_name. Using default time series: $selected_name."
        end

        # Collect all time series for the active device by looping through all weather years
        for (wy, data) in wind_ts_dict
            # Generate a unique time series name
            ts_name = "max_active_power_$(wy)"  # Unique name for each time series

            if selected_name in names(data)  # Check if the selected column exists in the DataFrame

                # Create the time series by extracting data from the DataFrame
                ts = SingleTimeSeries(;
                    name = ts_name,
                    data = TimeArray(tstamps, data[!, selected_name]), #create a TimeArray object w/ p.u. values
                    scaling_factor_multiplier = get_max_active_power, #apply scaling_factor_multiplier (we need this for SiennaPRASInterface runs)
                )

                # Store the time series along with its device name
                push!(ts_container, (device_name, ts))
            else
                @warn "Selected name $selected_name not found in data for year $year."
            end
        end
    end

    return ts_container
end

function add_wind_time_series_to_system!(sys::System, w_ts_container::Vector{Tuple{String, SingleTimeSeries}})
    # Loop through each wind RenewableDispatch object and add the time series
    for (device_name, time_series) in w_ts_container
        # Retrieve the active device by its name
        active_device = get_component(RenewableDispatch, sys, device_name)

        if active_device !== nothing
            # Add the time series to the system
            add_time_series!(sys, active_device, time_series)
            # println("Added time series: ", time_series.name, " to device: ", device_name)
        else
            @warn "Device $device_name not found in the system. Time series $(time_series.name) not added."
        end
    end
end

"""
note: not using create_split_load_forecasts function because inputs csvs are loading already represent the allocated load
therefore calling upload_split_load_forecasts function directly - see below
"""

#= function create_split_load_forecasts(load_dir::String, loads::Vector{PowerLoad})
    #################
    # Step 1: Read and Process Stochastic Loads
    #################
    # Read in the stochastic load values
    file_name = "Stochastic_Loads_CY2030_WY1998-2022.csv"
    file_path = joinpath(load_dir, file_name)
    df = CSV.read(file_path, DataFrame)

    # Filter out leap day (February 29)
    df = filter(row -> !(row.Month == 2 && row.Day == 29), df)

    # Initialize the dictionary to store each weather year as a separate DataFrame
    load_ts_dict = Dict{String, DataFrame}()

    # Define an 8760 vector based on a date interval object
    tstamps = collect(range(DateTime("2030-01-01T00:00:00"), DateTime("2030-12-31T23:00:00"), step = Dates.Hour(1)))

    # Populate the dictionary with weather year data
    for i in 1:25
        year_label = "Y$(1997 + i)"  # Construct year label (e.g., "Y1998", "Y1999", ...)
        load_values = df[!, Symbol(string(i))]  # Extract the load values for this year
        load_ts_dict[year_label] = DataFrame(:Time => tstamps, :Load => load_values)  # Add time and load columns
    end

    #################
    # Step 2: Read and Process Load Allocation Data
    #################
    # Read in monthly load allocations
    file_name = "NVE_Load_Share_Allocation.csv"
    file_path = joinpath(load_dir, file_name)
    df_demand_allocation = CSV.read(file_path, DataFrame)

    # Drop the :Month column
    df_demand_allocation = select!(df_demand_allocation, Not(:Month))

    #################
    # Step 3: Dynamically Initialize Split Load Dictionary
    #################
    # Initialize the dictionary to store split load time series
    split_load_ts_dict = Dict{String, Dict{String, DataFrame}}()

    # Populate split_load_ts_dict` with PowerLoad names
    for load in loads
        load_name = get_name(load)  # Dynamically retrieve the name of each PowerLoad
        split_load_ts_dict[load_name] = Dict{String, DataFrame}()  # Initialize an inner dictionary for each load
    end

    #################
    # Step 4: Populate Split Load Time Series
    #################
    for (wy, load_df) in load_ts_dict
        # Extract hourly load data and timestamps
        timestamps = load_df[:, :Time]
        total_load = load_df[:, :Load]

        # Iterate through each PowerLoad in df_demand_allocation
        for (col, load_name) in zip(eachcol(df_demand_allocation), names(df_demand_allocation))
            # Initialize an array for the split load
            split_load = Vector{Float64}(undef, length(total_load))

            # Split the load based on the allocation factors
            for i in 1:length(total_load)
                current_month = month(timestamps[i])  # Extract the current month (1–12)
                allocation_factor = col[current_month]  # Get the allocation factor for this month
                split_load[i] = total_load[i] * allocation_factor  # Apply the allocation factor
            end

            # Create a DataFrame for the split load
            split_load_df = DataFrame(:Time => timestamps, :Load => split_load)

            # Store the DataFrame in the appropriate dictionary
            split_load_ts_dict[load_name][wy] = split_load_df
        end
    end

    #################
    # Step 5: Return Results
    #################
    return load_ts_dict, split_load_ts_dict
end =#

function upload_split_load_forecasts(load_dir::String, loads::Vector{PowerLoad})
    #################
    # Step 1: Read and Process Stochastic regional loads into dataframes; store in separate dictionaries
    #################
    # Read in the stochastic load file (nevada power)
    np_power_file_name = "Stochastic_Loads_CY2030_WY1998-2022_Nevada_Power.csv"
    file_path_np = joinpath(load_dir, np_power_file_name)
    df_np = CSV.read(file_path_np, DataFrame)
    df_np = filter(row -> !(row.Month == 2 && row.Day == 29), df_np) # Filter out leap day (February 29)

    # Read in the stochastic load file (sierra)
    sierra_power_file_name = "Stochastic_Loads_CY2030_WY1998-2022_Sierra.csv"
    file_path_sierra = joinpath(load_dir, sierra_power_file_name)
    df_sierra = CSV.read(file_path_sierra, DataFrame)
    df_sierra = filter(row -> !(row.Month == 2 && row.Day == 29), df_sierra) # Filter out leap day (February 29)

    #################
    # Step 2: Initialize Dictionaries for Each Weather Year
    #################
    # Initialize the dictionaries to store each weather year as a separate DataFrame
    load_ts_dict = Dict{String, DataFrame}()
    # Define an 8760 vector based on a date interval object
    tstamps = collect(range(DateTime("2030-01-01T00:00:00"), DateTime("2030-12-31T23:00:00"), step = Dates.Hour(1)));

    # Populate combined dictionary for both regions
    for i in 1:25 # Populate the dictionaries with weather year data (Y1998 -> 1; Y2022 -> 25)
        year_label = "Y$(1997 + i)"  # Construct year label (e.g., "Y1998", "Y1999", ...)
        load_values_np = df_np[!, Symbol(string(i))]  # Nevada Power load values
        load_values_sierra = df_sierra[!, Symbol(string(i))]  # Sierra load values
        load_ts_dict["Nevada Power $year_label"] = DataFrame(:Time => tstamps, :Load => load_values_np)
        load_ts_dict["Sierra $year_label"] = DataFrame(:Time => tstamps, :Load => load_values_sierra)
    end

    load_ts_dict

    #################
    # Step 3: Convert DataFrames to PowerSystems.jl TimeSeries
    #################
    # initialize container for storing timeseries
    ts_container = Vector{Tuple{String, SingleTimeSeries}}()

    for load_object in loads

        # testing / troubleshooting
        # load_object = loads[1]

        # get device name
        device_name = get_name(load_object)

        # get max_active_power
        max_load = get_max_active_power(load_object)

        # Filter relevant time series for the device name
        relevant_ts = filter(k -> startswith(k, device_name), keys(load_ts_dict))

        for forecast_label in relevant_ts

            # testing / troubleshooting
            # forecast_label = first(relevant_ts)

            year_component = split(forecast_label)[end]  #split the string and take just the year component
            ts_name = "max_active_power_$(year_component)"
            normalized_data = load_ts_dict[forecast_label][:, :Load] ./ max_load

            ts = SingleTimeSeries(;
                name = ts_name,
                data = TimeArray(tstamps, normalized_data), # apply p.u. load forecast_label
                scaling_factor_multiplier = get_max_active_power,
            )

            # Store the time series and device name in the combined container
            push!(ts_container, (device_name, ts))
        end
    end

    #################
    # Step 4: Return Results
    #################
    return ts_container
end

function add_load_time_series_to_system!(sys::System, load_fx_timeseries::Vector{Tuple{String, SingleTimeSeries}})
    for (device_name, time_series) in load_fx_timeseries
        # Retrieve the active device by its name
        active_device = get_component(PowerLoad, sys, device_name)

        if active_device !== nothing
            # Add the time series to the system
            add_time_series!(sys, active_device, time_series)
            println("Added time series: ", time_series.name, " to device: ", device_name)
        else
            @warn "Device $device_name not found in the system. Time series $(time_series.name) not added."
        end
    end
end

#Function to update max_active_power time series for loads for SiennaPRASInterface runs
function update_max_active_power_ts!(sys, load_name, ts_name)
    # Retrieve the active load by its name
    active_load = get_component(PowerLoad, sys, load_name)
    #retrieve time series data for weather year of interest
    ts = get_time_series_array(SingleTimeSeries, active_load, ts_name; ignore_scaling_factors = true)
    #create new time series object using same name as deterministic time series
    ts_d = SingleTimeSeries(name = "max_active_power", data = ts)
    # remove existing / old time series definition for max active power from the system 
    remove_time_series!(sys, SingleTimeSeries, active_load, "max_active_power")
    # add new time series with updated definition for max active power to the system
    add_time_series!(sys, active_load, ts_d)
end
