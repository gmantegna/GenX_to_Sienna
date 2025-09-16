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

    data_root_dir = "/scratch/gpfs/gm1710/Sienna_cases"

    # Define key directories
    data_dir = joinpath(data_root_dir, "Sonoma", "GENX_Output") # GenX input files
    PSI_results_dir = joinpath(data_root_dir, "Sonoma", "Sienna_Outputs","output_results","PowerSimulations") # PSI processed results
    PRAS_results_dir = joinpath(data_root_dir, "Sonoma", "Sienna_Outputs","output_results","PRAS") # PRAS processed results
    sienna_simulation_dir = joinpath(data_root_dir, "Sonoma", "Sienna_Outputs","simulation_files") # Sienna simulation results
    
    # Ensure PSI output directory exists
    if !ispath(PSI_results_dir)
        mkpath(PSI_results_dir)
    end

    # Ensure PRAS output directory exists
    if !ispath(PRAS_results_dir)
        mkpath(PRAS_results_dir)
    end

    # Ensure Sienna simulation directory exists
    if !ispath(sienna_simulation_dir)
        mkpath(sienna_simulation_dir)
    end

    # Return all paths as a dictionary for easy access
    return Dict(
        :root_dir => root_dir,
        :data_dir => data_dir,
        :PSI_results_dir => PSI_results_dir,
        :PRAS_results_dir => PRAS_results_dir,
        :sienna_simulation_dir => sienna_simulation_dir
    )
end

function process_demand_data(demand_data_path::String, zone_dict::OrderedDict{String,Int64})
    # Read in the original GENX demand data
    demand_df = CSV.read(demand_data_path, DataFrame)

    # Create a mapping from zone numbers to zone names
    zone_number_to_name = Dict{Int,String}()
    for (name, number) in zone_dict
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

    # First extend the DataFrame by duplicating the last 24 entries
    last_24_rows = demand_df[end-23:end, :]
    demand_df = vcat(demand_df, last_24_rows)
    
    # Generate datetime range and add to DataFrame
    start_date = DateTime(1998, 1, 1)
    end_date = DateTime(2020, 12, 31, 23)
    dates = collect(start_date:Dates.Hour(1):end_date)
    
    # Add datetime and index columns to DataFrame
    demand_df.DateTime = dates
    demand_df.Y_index = Dates.year.(dates)
    demand_df.M_index = Dates.month.(dates)
    demand_df.D_index = Dates.day.(dates)
    
    #rearrange column order
    demand_df = select!(demand_df, :Y_index, :M_index, :D_index, :DateTime, Not([:Y_index, :M_index, :D_index, :DateTime]))
    
    # Filter out leap days (February 29th) using subset
    demand_df = subset(demand_df, [:M_index, :D_index] => (x, y) -> .!((x .== 2) .& (y .== 29)))

    return demand_df
end

function process_fuel_data(fuels_df::DataFrame)
    # Drop the first row (CO2 emissions) and the first column (index)
    fuels_df = fuels_df[2:end, 2:end]
    
    # First extend the DataFrame by duplicating the last 24 entries
    last_24_rows = fuels_df[end-23:end, :]
    fuels_df = vcat(fuels_df, last_24_rows)
    
    # Generate datetime range and add to DataFrame
    start_date = DateTime(1998, 1, 1)
    end_date = DateTime(2020, 12, 31, 23)
    dates = collect(start_date:Dates.Hour(1):end_date)
    
    # Add datetime and index columns to DataFrame
    fuels_df.DateTime = dates
    fuels_df.Y_index = Dates.year.(dates)
    fuels_df.M_index = Dates.month.(dates)
    fuels_df.D_index = Dates.day.(dates)
    
    # Rearrange column order
    fuels_df = select!(fuels_df, :Y_index, :M_index, :D_index, :DateTime, Not([:Y_index, :M_index, :D_index, :DateTime]))
    
    # Filter out leap days (February 29th) using subset
    fuels_df = subset(fuels_df, [:M_index, :D_index] => (x, y) -> .!((x .== 2) .& (y .== 29)))

    return fuels_df
end

function process_generator_variability_data(gen_var_data_path::String)
    # Read the generator variability data
    gen_var_df = CSV.read(gen_var_data_path, DataFrame)
    
    # Drop the first column (time index)
    select!(gen_var_df, Not(1))
    
    # First extend the DataFrame by duplicating the last 24 entries
    last_24_rows = gen_var_df[end-23:end, :]
    gen_var_df = vcat(gen_var_df, last_24_rows)
    
    # Generate datetime range and add to DataFrame
    start_date = DateTime(1998, 1, 1)
    end_date = DateTime(2020, 12, 31, 23)
    dates = collect(start_date:Dates.Hour(1):end_date)
    
    # Add datetime and index columns to DataFrame
    gen_var_df.DateTime = dates
    gen_var_df.Y_index = Dates.year.(dates)
    gen_var_df.M_index = Dates.month.(dates)
    gen_var_df.D_index = Dates.day.(dates)
    
    # Rearrange column order
    gen_var_df = select!(gen_var_df, :Y_index, :M_index, :D_index, :DateTime, Not([:Y_index, :M_index, :D_index, :DateTime]))
    
    # Filter out leap days (February 29th) using subset
    gen_var_df = subset(gen_var_df, [:M_index, :D_index] => (x, y) -> .!((x .== 2) .& (y .== 29)))

    return gen_var_df
end

function process_hydro_budget_data(hydro_budget_data_path::String)
    # Read the hydro budget data
    hydro_budget_df = CSV.read(hydro_budget_data_path, DataFrame)

    # drop first column (time index)
    select!(hydro_budget_df, Not(1))
    
    # First extend the DataFrame by duplicating the last 24 entries
    last_24_rows = hydro_budget_df[end-23:end, :]
    hydro_budget_df = vcat(hydro_budget_df, last_24_rows)
    
    # Generate datetime range and add to DataFrame
    start_date = DateTime(1998, 1, 1)
    end_date = DateTime(2020, 12, 31, 23)
    dates = collect(start_date:Dates.Hour(1):end_date)
    
    # Add datetime and index columns to DataFrame
    hydro_budget_df.DateTime = dates
    hydro_budget_df.Y_index = Dates.year.(dates)
    hydro_budget_df.M_index = Dates.month.(dates)
    hydro_budget_df.D_index = Dates.day.(dates)
    
    # Rearrange column order
    hydro_budget_df = select!(hydro_budget_df, :Y_index, :M_index, :D_index, :DateTime, Not([:Y_index, :M_index, :D_index, :DateTime]))
    
    # Filter out leap days (February 29th) using subset
    hydro_budget_df = subset(hydro_budget_df, [:M_index, :D_index] => (x, y) -> .!((x .== 2) .& (y .== 29)))

    return hydro_budget_df
end

function create_storage_parameters_df(Storage_objects::Vector{EnergyReservoirStorage}, output_path::String)
    # Create DataFrame with column names matching the parameters we want to check
    df_storage = DataFrame(
        name = String[],
        base_power = Float64[],
        storage_capacity = Float64[],
        storage_level_limits_min = Float64[],
        storage_level_limits_max = Float64[],
        initial_storage_level = Float64[],
        efficiency_in = Float64[],
        efficiency_out = Float64[],
        input_power_limits_min = Float64[],
        input_power_limits_max = Float64[],
        output_power_limits_min = Float64[],
        output_power_limits_max = Float64[],
        rating = Float64[]
    )

    # Loop through all storage devices and add their parameters
    for storage in Storage_objects
        push!(df_storage, (
            get_name(storage),
            get_base_power(storage),
            get_storage_capacity(storage),
            get_storage_level_limits(storage).min,
            get_storage_level_limits(storage).max,
            get_initial_storage_capacity_level(storage),
            get_efficiency(storage).in,
            get_efficiency(storage).out,
            get_input_active_power_limits(storage).min,
            get_input_active_power_limits(storage).max,
            get_output_active_power_limits(storage).min,
            get_output_active_power_limits(storage).max,
            get_rating(storage)
        ))
    end

    # Write the DataFrame to a CSV file
    CSV.write(output_path, df_storage)

end

function create_area_interchange_parameters_df(area_interchanges::Vector{AreaInterchange}, output_path::String)
    # Create DataFrame with column names matching the parameters we want to check
    df_interchanges = DataFrame(
        name = String[],
        from_area = String[],
        to_area = String[],
        flow_limits_min = Float64[],
        flow_limits_max = Float64[]
    )

    # Loop through all area interchanges and add their parameters
    for interchange in area_interchanges
        flow_limits = get_flow_limits(interchange)
        push!(df_interchanges, (
            get_name(interchange),
            get_name(get_from_area(interchange)),
            get_name(get_to_area(interchange)),
            get_flow_limits(interchange).to_from,
            get_flow_limits(interchange).from_to
        ))
    end

    # Write the DataFrame to a CSV file
    CSV.write(output_path, df_interchanges)
end

function create_line_parameters_df(lines::Vector{Line}, output_path::String)
    # Create DataFrame with column names matching the parameters we want to check
    df_lines = DataFrame(
        name = String[],
        from_bus = String[],
        to_bus = String[],
        rating = Float64[]
    )

    # Loop through all lines and add their parameters
    for line in lines
        arc = get_arc(line)
        push!(df_lines, (
            get_name(line),
            get_name(get_from(arc)),
            get_name(get_to(arc)),
            get_rating(line)
        ))
    end

    # Write the DataFrame to a CSV file
    CSV.write(output_path, df_lines)
end

function create_powerload_parameters_df(sys::System, paths::Dict)
    # Initialize DataFrame with columns for PowerLoad attributes
    df = DataFrame(
        name = String[],
        base_power = Float64[],
        active_power = Float64[],
        bus_name = String[]
    )

    # Loop through all PowerLoad components in the system
    for load in get_components(PowerLoad, sys)
        # Get the attributes
        name = get_name(load)
        base_power = get_base_power(load)
        active_power = get_active_power(load)
        bus = get_bus(load)
        bus_name = get_name(bus)

        # Add row to DataFrame
        push!(df, (name, base_power, active_power, bus_name))
    end

    # Write to CSV
    CSV.write(joinpath(paths[:data_dir], "powerload_parameters.csv"), df)
end

function create_transmission_interface_parameters_df(sys::System, paths::Dict)
    # Initialize DataFrame with columns for TransmissionInterface attributes
    df = DataFrame(
        name = String[],
        active_power_flow_limits_min = Float64[],
        active_power_flow_limits_max = Float64[],
        direction_mapping = Dict{String, Int}[]
    )

    # Loop through all TransmissionInterface components in the system
    for interface in get_components(TransmissionInterface, sys)
        # Get the attributes
        name = get_name(interface)
        flow_limits = get_active_power_flow_limits(interface)
        direction_mapping = get_direction_mapping(interface)

        # Add row to DataFrame
        push!(df, (name, flow_limits.min, flow_limits.max, direction_mapping))
    end

    # Write to CSV
    CSV.write(joinpath(paths[:data_dir], "transmission_interface_parameters.csv"), df)
end

#= function create_PHS_parameters_df(pumped_hydro_objects::Vector{HydroPumpedStorage}, output_path::String)
    # Create DataFrame with column names matching the parameters we want to check
    df_pumped_hydro = DataFrame(
        name = String[],
        base_power = Float64[],
        storage_capacity_up = Float64[],
        storage_capacity_down = Float64[],
        initial_storage_up = Float64[],
        initial_storage_down = Float64[],
        pump_efficiency = Float64[],
        active_power_limits_min = Float64[],
        active_power_limits_max = Float64[],
        rating = Float64[]
    )

    # Loop through all pumped hydro devices and add their parameters
    for ph in pumped_hydro_objects
        push!(df_pumped_hydro, (
            get_name(ph),
            get_base_power(ph),
            get_storage_capacity(ph).up,
            get_storage_capacity(ph).down,
            get_initial_storage(ph).up,
            get_initial_storage(ph).down,
            get_pump_efficiency(ph),
            get_active_power_limits(ph).min,
            get_active_power_limits(ph).max,
            get_rating(ph)
        ))
    end

    # Write the DataFrame to a CSV file
    CSV.write(output_path, df_pumped_hydro)
    
    return df_pumped_hydro
end =#

function create_PHS_storage_parameters_df(PHS_Storage_objects::Vector{EnergyReservoirStorage}, output_path::String)
    # Create DataFrame with column names matching the parameters we want to check
    df_PHS_storage = DataFrame(
        name = String[],
        base_power = Float64[],
        storage_capacity = Float64[],
        storage_level_limits_min = Float64[],
        storage_level_limits_max = Float64[],
        initial_storage_level = Float64[],
        efficiency_in = Float64[],
        efficiency_out = Float64[],
        input_power_limits_min = Float64[],
        input_power_limits_max = Float64[],
        output_power_limits_min = Float64[],
        output_power_limits_max = Float64[],
        rating = Float64[]
    )

    # Loop through all storage devices and add their parameters
    for PHS_storage in PHS_Storage_objects
        push!(df_PHS_storage, (
            get_name(PHS_storage),
            get_base_power(PHS_storage),
            get_storage_capacity(PHS_storage),
            get_storage_level_limits(PHS_storage).min,
            get_storage_level_limits(PHS_storage).max,
            get_initial_storage_capacity_level(PHS_storage),
            get_efficiency(PHS_storage).in,
            get_efficiency(PHS_storage).out,
            get_input_active_power_limits(PHS_storage).min,
            get_input_active_power_limits(PHS_storage).max,
            get_output_active_power_limits(PHS_storage).min,
            get_output_active_power_limits(PHS_storage).max,
            get_rating(PHS_storage)
        ))
    end

    # Write the DataFrame to a CSV file
    CSV.write(output_path, df_PHS_storage)

end

function query_write_export_results(sim::Simulation, path_scenario::String, uc_decision_name::String)
    ###########################
    # Query Results
    ###########################
    sim_results = SimulationResults(sim)
    results = get_decision_problem_results(sim_results, uc_decision_name) # UC stage result metadata

    # Input TimeSeries Parameters
    load_parameter = read_realized_parameter(results, "ActivePowerTimeSeriesParameter__PowerLoad")
    thermal_parameter = read_realized_parameter(results, "ActivePowerTimeSeriesParameter__ThermalStandard")
    renewDispatch_parameter = read_realized_parameter(results, "ActivePowerTimeSeriesParameter__RenewableDispatch")
    renewNonDispatch_parameter = read_realized_parameter(results, "ActivePowerTimeSeriesParameter__RenewableNonDispatch") # Sienna doesnt store nonDispatch
    hydro_parameter = read_realized_parameter(results, "ActivePowerTimeSeriesParameter__HydroDispatch")

    # Output Realized Generation Values
    thermal_active_power = read_realized_variable(results, "ActivePowerVariable__ThermalStandard")
    renewDispatch_active_power = read_realized_variable(results, "ActivePowerVariable__RenewableDispatch")
    # renewNonDispatch_active_power = read_realized_variable(results, "ActivePowerVariable__RenewableNonDispatch") Q: why isnt this available in results; bc no dispatch?
    hydro_active_power = read_realized_variable(results, "ActivePowerVariable__HydroDispatch")
    storage_charge = read_realized_variable(results, "ActivePowerInVariable__EnergyReservoirStorage")
    storage_discharge = read_realized_variable(results, "ActivePowerOutVariable__EnergyReservoirStorage")

    # combine all FTM generators
    gen_active_power = hcat(thermal_active_power, select(renewDispatch_active_power, Not(1)), select(hydro_active_power, Not(1)))

    # Output Realized TX flows
    tx_flow = read_realized_variable(results, "FlowActivePowerVariable__AreaInterchange")

    # Output Expressions
    power_balance = read_realized_expression(results, "ActivePowerBalance__Area")

    # Get Production Costs
    pc_thermal = read_realized_expression(results, "ProductionCostExpression__ThermalStandard")
    pc_renewable = read_realized_expression(results, "ProductionCostExpression__RenewableDispatch")
    pc_hydro = read_realized_expression(results, "ProductionCostExpression__HydroDispatch")
    #all_pc = hcat(pc_thermal,select(pc_renewable, Not(1)))
    all_pc = hcat(pc_thermal,select(pc_renewable, Not(1)), select(pc_hydro, Not(1)))

    ###########################
    # Export Results
    ###########################
    # Define output paths and write dataframes to CSV
    CSV.write(joinpath(path_scenario, "load_active_power.csv"), load_parameter) # Input time series values
    CSV.write(joinpath(path_scenario, "FTM_renewable_parameters.csv"), renewDispatch_parameter) # Input time series values
    CSV.write(joinpath(path_scenario, "thermal_parameters.csv"), thermal_parameter) # Input time series values
    CSV.write(joinpath(path_scenario, "hydro_parameter.csv"), hydro_parameter) # Input time series values
    CSV.write(joinpath(path_scenario, "BTM_active_power.csv"), renewNonDispatch_parameter) # Input time series values

    CSV.write(joinpath(path_scenario, "FTM_renewable_active_power.csv"), renewDispatch_active_power)
    CSV.write(joinpath(path_scenario, "thermal_active_power.csv"), thermal_active_power)
    CSV.write(joinpath(path_scenario, "FTM_generator_active_power.csv"), gen_active_power)
    CSV.write(joinpath(path_scenario, "tx_flow.csv"), tx_flow)
    CSV.write(joinpath(path_scenario, "storage_charge.csv"), storage_charge)
    CSV.write(joinpath(path_scenario, "storage_discharge.csv"), storage_discharge)
    CSV.write(joinpath(path_scenario, "power_balance.csv"), power_balance)
    CSV.write(joinpath(path_scenario, "production_costs.csv"), all_pc)    
end

function system_capacity_query(unit_collection::Dict, paths::Dict)
    # Initialize an empty DataFrame
    df = DataFrame(Resource = String[], MW_capacity = Float64[])

    # for each generator collection, loop through each unit w/in that collection
    for (category, gen_collection) in unit_collection
        if category == "StorageUnits" #batteries
            for unit in gen_collection
                if get_available(unit)  # Check if the unit is active
                    name = get_name(unit)  # Get the generator's name
                    capacity = get_output_active_power_limits(unit).max  # Get the max active discharge power (MW)

                    # Append to DataFrame
                    push!(df, (name, capacity))
                else
                # do nothing
                end
            end
        else # all other generator types
            for unit in gen_collection
                if get_available(unit)  # Check if the unit is active

                    name = get_name(unit)  # Get the generator's name
                    capacity = get_max_active_power(unit)  # Get the max active power (MW)

                    # Append to DataFrame
                    push!(df, (name, capacity))
                else
                    #do nothing
                end
            end
        end # if loop
    end

    #write df to csv
    CSV.write(joinpath(paths[:data_dir], "nameplate_capacity.csv"), df);
end
