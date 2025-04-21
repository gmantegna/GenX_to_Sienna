function initialize_paths_and_inputs(r2x_output_name::String, scenario_name::String)
    # Define root directory and related paths
    current_dir = pwd()

    # Check if "Projects" exists in the path (want to have root directory get you to \\FPA_Sienna)
    if occursin("Projects", current_dir)
        # Find the position where "Projects" starts
        start_pos = first(findfirst("Projects", current_dir)) #return the starting index of the string position for "Projects"
        # Trim off everything from "Projects" onward
        root_dir = current_dir[1:start_pos-1]
    else
    end
    # @show root_dir

    data_dir = joinpath(root_dir, "Projects", "NVE") # root data directory
    r2x_dir = joinpath(data_dir, r2x_output_name) # R2X folder
    output_dir = joinpath(root_dir, "Projects", "NVE", "sienna_runs", "run_output") # Sienna processed outputs and simulation files
    sim_files_dir = joinpath(output_dir, "simulation_files") # Sienna outputs

    # Directory for output results
    scenario_dir_deterministic = joinpath(output_dir, scenario_name,"deterministic") # Scenario name passed as argument
    scenario_dir_stochastics = joinpath(output_dir, scenario_name,"stochastics") # Scenario name passed as argument

    # Stochastic inputs
    FTM_solar_dir = joinpath(data_dir, "LOLP_inputs", "solar") #FTM
    BTM_solar_dir = joinpath(data_dir, "LOLP_inputs", "DPV") #BTM
    wind_dir = joinpath(data_dir, "LOLP_inputs", "wind")
    load_dir = joinpath(data_dir, "LOLP_inputs", "load")
    LOLP_inputs = joinpath(data_dir, "LOLP_inputs")

    # Ensure that simulation_files exists
    if !ispath(sim_files_dir)
        mkpath(sim_files_dir)
    else
    end

    # Ensure that deterministic directory exists
    if !ispath(scenario_dir_deterministic)
        mkpath(scenario_dir_deterministic)
    else
    end

    # Ensure that stochastic directory exists
    if !ispath(scenario_dir_stochastics)
        mkpath(scenario_dir_stochastics)
    else
    end

    # Point to the appropriate input file
    descriptors_dir = joinpath(root_dir, "Projects", "NVE", "sienna_runs", "user_descriptors")

    # Return all paths as a dictionary for easy access
    return Dict(
        :root_dir => root_dir,
        :data_dir => data_dir,
        :r2x_dir => r2x_dir,
        :output_dir => output_dir,
        :sim_files_dir => sim_files_dir,
        :scenario_dir_d => scenario_dir_deterministic,
        :scenario_dir_s => scenario_dir_stochastics,
        :FTM_solar_dir => FTM_solar_dir,
        :BTM_solar_dir => BTM_solar_dir,
        :wind_dir => wind_dir,
        :load_dir => load_dir,
        :LOLP_inputs => LOLP_inputs,
        :descriptors_dir => descriptors_dir
    )
end

function build_system_from_r2x(paths::Dict{Symbol, String}, r2x_output_name::String)
    # Base power for the system
    base_power = 1.0

    # Extract paths from the dictionary
    root_dir = paths[:root_dir]
    # data_dir = paths[:data_dir]

    # Define paths to input files
    descriptors = joinpath(paths[:descriptors_dir], "$r2x_output_name.yaml")
    generator_mapping = joinpath(root_dir, "Projects", "NVE", "sienna_runs", "generator_mapping.yaml")
    timeseries_metadata = joinpath(root_dir, "Projects", "NVE", r2x_output_name, "timeseries_pointers.json")

    # Load the data
    data = PowerSystemTableData(
        paths[:r2x_dir],
        base_power,
        descriptors;
        generator_mapping_file = generator_mapping,
        timeseries_metadata_file = timeseries_metadata,
    )

    #call the helper function modify_data
    modify_data(data)

    # Create the system
    sys = System(data, time_series_in_memory= true) #Q: using from memory (smaller systems; faster result) of hardrive sys what is time_series_in_memory?

    # Define the path to the PWL output CSV
    pw_file = joinpath(paths[:r2x_dir], "generator_PWL_output.csv")

    # Read the PWL output into a DataFrame
    pw_data = DataFrame(CSV.File(pw_file))

    # Return the system and the PWL data
    return sys, pw_data
end

function modify_system!(sys::System, pw_data::DataFrame, paths::Dict)
    #create area interchanges so that system goes from 1 area to 4
    create_area_interchanges(sys, paths[:r2x_dir])
    # attach_reserve_requirements(sys, joinpath(paths[data_dir], r2x_output_name))

    # Define initial state of charge (SOC) for storage
    for storage in get_components(EnergyReservoirStorage, sys)
        set_initial_storage_capacity_level!(storage, 0.33)
    end

    # Set generator services (i.e., assign memberships)
    for gen in get_components(ThermalStandard, sys)
        if !occursin("Purchases", get_name(gen))
            #concat the two vectors of services
            eligible_services = vcat(
                collect(get_components(VariableReserve{ReserveUp}, sys)),
                # disabling non-spinning reserves Q: what is going on here
                # collect(get_components(VariableReserveNonSpinning, sys))
                )
            set_services!(gen, eligible_services)
        end
    end

    # Disable all reserves
    for reserve in get_components(VariableReserve{ReserveUp}, sys)
        set_available!(reserve, false)
    end

    # Constrain the Market imports (i.e., convert mkt purchases to renewable dispatch)
    #constrain_market_imports(sys, paths[:r2x_dir])

    # Set market purchases availability
    # set_available!(get_component(ThermalStandard, sys, "Southern Purchases (NVP)"), false)
    # set_available!(get_component(ThermalStandard, sys, "Northern Purchases (Sierra)"), false)
    # set_available!(get_component(RenewableDispatch, sys, "Sierra Solar II"), false) # not active in PLEXOS system
    # set_available!(get_component(EnergyReservoirStorage, sys, "Sierra Solar II BESS"), false) # not active in PLEXOS system

    # Modify generator heat rates
    for row in eachrow(pw_data)
        generator_name = row.name
        generator = get_component(ThermalStandard, sys, generator_name)
        generator_op_cost = get_operation_cost(generator)

        new_vc = PiecewisePointCurve(
            [
                (row.output_point_0 , row.heat_rate_avg_0),
                (row.output_point_1 , row.heat_rate_incr_1),
                (row.output_point_2 , row.heat_rate_incr_2),
                (row.output_point_3 , row.heat_rate_incr_3),
                (row.output_point_4 , row.heat_rate_incr_4),
                (row.output_point_5 , row.heat_rate_incr_5),
            ]
        )

        fuel_curve = FuelCurve(
            value_curve = new_vc,
            fuel_cost = generator_op_cost.variable.fuel_cost,
            vom_cost = generator_op_cost.variable.vom_cost
            )

        new_generator_cost = ThermalGenerationCost(
            variable = fuel_curve,
            fixed = 0,
            start_up = generator_op_cost.start_up,
            shut_down = generator_op_cost.shut_down
            )
        set_operation_cost!(generator, new_generator_cost)
    end

    #reassign ThermalStandard units that get flagged as NG
    geothermal_units = ["Gerlach", "Lone Mountain", "North Valley 2", "North Valley Eavor Loop", "Pinto"]
    for unit in geothermal_units
        #check to see if fuel type is geothermal
        active_geo = get_component(ThermalStandard, sys, unit)
        if get_fuel(active_geo) == ThermalFuels.GEOTHERMAL
            # do nothing
        else
            set_fuel!(active_geo, ThermalFuels.GEOTHERMAL)
            #set the fuel cost to 0
            set_fuel_cost!(sys, active_geo, 0.0)
        end
    end

#=     for unit in get_components(ThermalStandard, sys)
        if get_fuel(unit) == ThermalFuels.GEOTHERMAL
            println("Geothermal Unit: ", get_name(unit))
        end
    end =#

    # Set the reference (slack) bus
    set_bustype!(get_component(ACBus, sys, "Nevada Power"), "REF");
end

function Tuning_Adjustment_Costs!(sys::System)
    # review Tuning Adjustment operational cost
    active_unit = get_component(ThermalStandard, sys, "Tuning Adjustment");
    get_operation_cost(active_unit);

    # set availability
    # set_available!(get_component(ThermalStandard, sys, "Tuning Adjustment"), false)

    # define $2000/MWh cost CostCurve (w/ no VOM)
    TA_CC = CostCurve(LinearCurve(2000),UnitSystem.NATURAL_UNITS,LinearCurve(0));
    # define Operational Cost via ThermalGenerationCost template
    op_cost = ThermalGenerationCost(; variable = TA_CC, fixed = 0, start_up = 0, shut_down = 0,);
    # assign the operational cost to the Tuning Adjustment unit
    set_operation_cost!(get_component(ThermalStandard, sys, "Tuning Adjustment"), op_cost);
    # check your work
    get_operation_cost(active_unit);
end

function Demand_Response_CleanUp!(sys::System)

    # deactivate original DR units (R2X picks these up as curtailable resources bc of the "*PV" wildcard)
    set_available!(get_component(RenewableDispatch, sys, "NVE.Owned.DR/DSM.PV"), false)
    set_available!(get_component(RenewableDispatch, sys, "PPA.Contracted.DR/DSM.PV"), false)

    # Remove original DR units from the system as a precautionary to make sure PRAS isnt picking these up)
    remove_component!(RenewableDispatch, sys, "NVE.Owned.DR/DSM.PV")
    remove_component!(RenewableDispatch, sys, "PPA.Contracted.DR/DSM.PV")

    # define $3000/MWh cost CostCurve (w/ no VOM); clears after tuning adjustment per benchmark system
    DR_CC = CostCurve(LinearCurve(3000),UnitSystem.NATURAL_UNITS,LinearCurve(0))

    # instantiate new DR units
    dr_names = ["NVE_Owned_DR_DSM_ACLM", "NVE_Owned_DR_DSM_BESS-2", "NVE_Owned_DR_DSM_PV", "PPA_Contracted_DR_DSM_PV"]
    dr_buses = ["Nevada Power", "Nevada Power", "Nevada Power", "Nevada Power"]
    dr_ratings = [283, 4, 278, 4] # MW

    # collection for DR units we create
    dr_units = Dict()

    # Loop through the resources and create each ThermalStandard unit
    for (i, name) in enumerate(dr_names)
        dr_units[name] = ThermalStandard(
            name = name,
            available = true,
            status = true,
            bus = get_component(ACBus, sys, dr_buses[i]),
            active_power = 0.0, # Per-unitized by device base_power
            reactive_power = 0.0, # Per-unitized by device base_power
            rating = dr_ratings[i], # Per-unitized by device base_power
            active_power_limits = (min = 0.0, max = dr_ratings[i]), # Per-unitized by device base_power
            reactive_power_limits = nothing, # Per-unitized by device base_power
            ramp_limits = nothing, # we arent enforcing ramp limits per incumbent setup
            operation_cost = ThermalGenerationCost(; variable = DR_CC, fixed = 0, start_up = 0, shut_down = 0,),
            base_power = 1,
            time_limits = nothing, #we arent enforcing MDT/MUT per incumbent setup
            must_run = false,
            prime_mover_type = PrimeMovers.OT, # Other
            fuel = ThermalFuels.OTHER,
        )
    end

    # spot check a unit
    dr_units["NVE_Owned_DR_DSM_ACLM"]

    # add the units to the system
    for (name, unit) in dr_units
        add_component!(sys, unit)
        # println("Added ThermalStandard unit: $name")
    end

    # check your work
    active_DR_unit = get_component(ThermalStandard, sys, "NVE_Owned_DR_DSM_ACLM")
    # active_DR_unit = get_component(ThermalStandard, sys, "NVE_Owned_DR_DSM_BESS-2")
    get_operation_cost(active_DR_unit)
end

function convert_BTM_units!(sys::System, paths::Dict)

#Deactive initial BTM solar units bc R2X picks them up as curtailable resources)
set_available!(get_component(RenewableDispatch, sys, "DPV_Sierra"), false)
set_available!(get_component(RenewableDispatch, sys, "DPV_Nevada Power"), false)

#remove these units from system as a precautionary to make sure PRAS isnt picking these up
remove_component!(RenewableDispatch, sys, "DPV_Sierra")
remove_component!(RenewableDispatch, sys, "DPV_Nevada Power")   

# instantiate DPV_Nevada Power BTM
DPV_Nevada_Power_ND = RenewableNonDispatch(;
name = "DPV_Nevada Power ND",
available = true,
bus = get_component(ACBus, sys, "Nevada Power"),
active_power = 0.0, # Per-unitized by device base_power
reactive_power = 0.0, # Per-unitized by device base_power
rating = 1372.6, # Per-unitized by device base_power
prime_mover_type = PrimeMovers.PVe,
power_factor = 1.0,
base_power = 1,
);

# instantiate DPV_Sierra BTM
DPV_Sierra_ND = RenewableNonDispatch(;
name = "DPV_Sierra ND",
available = true,
bus = get_component(ACBus, sys, "Sierra"),
active_power = 0.0, # Per-unitized by device base_power
reactive_power = 0.0, # Per-unitized by device base_power
rating = 98.5, # Per-unitized by device base_power
prime_mover_type = PrimeMovers.PVe,
power_factor = 1.0,
base_power = 1,
);

# add BTM solar resources (rev) to the system
add_components!(sys, [DPV_Nevada_Power_ND, DPV_Sierra_ND])

# define and assign time series to the revised BTM solar units
input_RF_csv = joinpath(paths[:r2x_dir],"Data", "RenewableDispatch_Rating Factor_._2030.csv")
input_map_csv = joinpath(paths[:r2x_dir],"Data", "RenewableDispatch_max_active_power_._2030.csv")
df_RF = CSV.read(input_RF_csv, DataFrame)
df_map = CSV.read(input_map_csv, DataFrame)

# Extract and store time series for each DPV resource

# define DateTime
tstamps = df_RF.DateTime  # Use DateTime column directly

#DPV_Nevada Power (rev)
# create PowerSystems.jl timeseries for Rating Factor
ts_rf = SingleTimeSeries(;
    name = "Rating Factor",
    data = TimeArray(tstamps, df_RF[!, "DPV_Nevada Power"]./maximum(df_RF[!, "DPV_Nevada Power"])), #scaling by max value
    scaling_factor_multiplier = get_max_active_power,
    )

# create PowerSystems.jl timeseries for max_active_power
ts_map = SingleTimeSeries(;
    name = "max_active_power",
    data = TimeArray(tstamps, df_map[!, "DPV_Nevada Power"])./get_max_active_power(DPV_Nevada_Power_ND),
    scaling_factor_multiplier = get_max_active_power, #scaling by max value
    )

# add both timeseries to each BTM object
add_time_series!(sys, DPV_Nevada_Power_ND, ts_rf);
add_time_series!(sys, DPV_Nevada_Power_ND, ts_map);

# remove_time_series!(sys,SingleTimeSeries, DPV_Nevada_Power_ND, "Rating Factor")
# remove_time_series!(sys,SingleTimeSeries, DPV_Nevada_Power_ND, "max_active_power")

#= # check your work
active_unit = get_component(RenewableDispatch, sys, "DPV_Nevada Power")
active_BTM_unit = get_component(RenewableNonDispatch, sys, "DPV_Nevada Power ND")
get_time_series_array(SingleTimeSeries, active_unit, "Rating Factor"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_BTM_unit, "Rating Factor"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_unit,  "Rating Factor"; ignore_scaling_factors = false)
get_time_series_array(SingleTimeSeries, active_BTM_unit,  "Rating Factor"; ignore_scaling_factors = false)

get_time_series_array(SingleTimeSeries, active_unit, "max_active_power"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_BTM_unit, "max_active_power"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_unit,  "max_active_power"; ignore_scaling_factors = false)
get_time_series_array(SingleTimeSeries, active_BTM_unit,  "max_active_power"; ignore_scaling_factors = false) =#

#DPV_Sierra (rev)
# create PowerSystems.jl timeseries for Rating Factor
ts_rf = SingleTimeSeries(;
    name = "Rating Factor",
    data = TimeArray(tstamps, df_RF[!, "DPV_Sierra"]./maximum(df_RF[!, "DPV_Sierra"])), #scaling by max value
    scaling_factor_multiplier = get_max_active_power,
    )

# create PowerSystems.jl timeseries for max_active_power
ts_map = SingleTimeSeries(;
    name = "max_active_power",
    data = TimeArray(tstamps, df_map[!, "DPV_Sierra"])./get_max_active_power(DPV_Sierra_ND),
    scaling_factor_multiplier = get_max_active_power, #scaling by max value
    )

# add both timeseries to each BTM object
add_time_series!(sys, DPV_Sierra_ND, ts_rf);
add_time_series!(sys, DPV_Sierra_ND, ts_map);

# remove_time_series!(sys,SingleTimeSeries, DPV_Nevada_Power_ND, "Rating Factor")
# remove_time_series!(sys,SingleTimeSeries, DPV_Nevada_Power_ND, "max_active_power")

#= # check your work
active_unit = get_component(RenewableDispatch, sys, "DPV_Sierra")
active_BTM_unit = get_component(RenewableNonDispatch, sys, "DPV_Sierra ND")
get_time_series_array(SingleTimeSeries, active_unit, "Rating Factor"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_BTM_unit, "Rating Factor"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_unit,  "Rating Factor"; ignore_scaling_factors = false)
get_time_series_array(SingleTimeSeries, active_BTM_unit,  "Rating Factor"; ignore_scaling_factors = false)

get_time_series_array(SingleTimeSeries, active_unit, "max_active_power"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_BTM_unit, "max_active_power"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_unit,  "max_active_power"; ignore_scaling_factors = false)
get_time_series_array(SingleTimeSeries, active_BTM_unit,  "max_active_power"; ignore_scaling_factors = false) =#
end

function fix_Hydro_Dispatch!(sys::System)
    # review Tuning Adjustment operational cost
    active_hydro_unit = active_unit = get_component(HydroDispatch, sys,"Hoover Dam (NV)")

    # see what time series is already in the system
    show_time_series(active_hydro_unit)

    # remove the existing time series
    remove_time_series!(sys,SingleTimeSeries, active_hydro_unit, "max_active_power")
    remove_time_series!(sys,SingleTimeSeries, active_hydro_unit, "Max Energy Month")
    remove_time_series!(sys,SingleTimeSeries, active_hydro_unit, "Min Energy Hour")

    #pull in fixed dispatch profile from benchmark deterministic run_output
    hydro_input_csv = joinpath(paths[:r2x_dir],"Data", "hydro_ts_deterministic.csv")
    df_hydro = CSV.read(hydro_input_csv, DataFrame)

    # define PowerSystems.jl time series
    ## DateTime
    tstamps = df_hydro.DateTime  # Use DateTime column directly

    ## Instantiate Fixed Hydro Dispatch Profile SingleTimeSeries
    sts_fixed_hydro = SingleTimeSeries(;
        name = "max_active_power",
        data = TimeArray(tstamps, df_hydro[!, "Hoover Dam (NV)"])./get_max_active_power(active_hydro_unit), #normalize by nameplate capacity
        scaling_factor_multiplier = get_max_active_power, #scaling by max active power
    )

    ## now add the time series to the HydroDispatch object
    add_time_series!(sys, active_hydro_unit, sts_fixed_hydro);

    # check your work
    show_time_series(active_hydro_unit)
    get_time_series_array(SingleTimeSeries, active_hydro_unit, "max_active_power"; ignore_scaling_factors = true)
    get_time_series_array(SingleTimeSeries, active_hydro_unit, "max_active_power"; ignore_scaling_factors = false)

end

function ThermalStandard_missing_ts!(sys,paths::Dict{Symbol, String})
    # Extract paths from the dictionary
    root_dir = paths[:r2x_dir]

    # Define paths to input files
    file_path = joinpath(root_dir, "missing_fuel_timeseries.csv")

    # read in datafile
    df = CSV.read(file_path, DataFrame)

    # generate and assign fuel price timeseries to ThermalStandard generators that are missing them (dual fuel units)
    DFU = ["Valmy CT 3", "Valmy CT 4", "Silverhawk 3", "Silverhawk 4"]
    for device_name in DFU # loop through dual fuel units
        # Retrieve the active device by its name
        active_device = get_component(ThermalStandard, sys, device_name)

        # create the time series
        ts = SingleTimeSeries(
            name = "fuel_price",
            data = TimeArray(df[!,"DateTime"],df[!,device_name]), #create a TimeArray object
        )

        # Add the time series to the system
        add_time_series!(sys, active_device, ts)

        @show "Added fuel price time series to $device_name"
    end
end

function update_thermal_fuel_price_timeseries!(sys)
    # look at an example unit
    gen = get_component(ThermalStandard, sys, "Chuck Lenzi 1_A")
    ts_test = get_time_series(DeterministicSingleTimeSeries, gen, "fuel_price")
    get_time_series_array(DeterministicSingleTimeSeries, gen, "fuel_price")

    # retrieve forecast horizon and interval
    horizon = get_horizon(ts_test)
    interval = get_interval(ts_test)

    # remove DeterministicSingleTimeSeries (i.e. forecasts) from ALL objects (ThermalStandard, RenewableDispatch, etc.)
    remove_time_series!(sys, DeterministicSingleTimeSeries)

    #reassign fuel_price timeseries to ThermalStandard objects
    for g in collect(get_components(ThermalStandard, sys))
        if "fuel_price" ∈ get_name.(get_time_series_keys(g))
            fuel_ts = get_time_series(SingleTimeSeries, g, "fuel_price")
            fuel_array = get_time_series_array(SingleTimeSeries, g, "fuel_price"; ignore_scaling_factors = true)
            tstamp = timestamp(fuel_array)
            vals = values(fuel_array)
            new_ts = SingleTimeSeries("fuel_price", TimeArray(tstamp, vals))
            remove_time_series!(sys, SingleTimeSeries, g, "fuel_price")
            set_fuel_cost!(sys, g, new_ts)
        end
    end
    # add back in the forecasts for all objects (e.g., ThermalStandard, RenewableDispatch, etc.)
    transform_single_time_series!(sys, horizon, interval)
end

function define_device_model_non_weather(template_uc)
    # Define custom DeviceModel for ThermalStandard with time-varying max_active_power
    thermal_model = DeviceModel(ThermalStandard, ThermalStandardUnitCommitment; time_series_names = Dict{Any, String}(
                    PowerSimulations.FuelCostParameter => "fuel_price",
                    PowerSimulations.ActivePowerTimeSeriesParameter => "max_active_power",))
    PowerSimulations.set_device_model!(template_uc, thermal_model)

    # Define Hydro model
    PowerSimulations.set_device_model!(template_uc, HydroDispatch, HydroDispatchRunOfRiver)

    # Define custom DeviceModel for EnergyReservoirStorage
    storage_model = DeviceModel(
        EnergyReservoirStorage,
        StorageDispatchWithReserves;
        attributes = Dict(
            "reservation" => true, # True prevents discharging and charging in the same period
            "cycling_limits" => false,
            "energy_target" => false,
            "complete_coverage" => false,
            "regularization" => true, # Regularizes storage dispatch to prevent large swings in dispatch.
        ),
    )
    # Assign the storage model to the template_uc
    PowerSimulations.set_device_model!(template_uc, storage_model)
end

function define_device_model_weather_deterministic(template_uc)
    # define the load model
    PowerSimulations.set_device_model!(template_uc, PowerLoad, StaticPowerLoad)

    # define the renewable dispatch model
    PowerSimulations.set_device_model!(template_uc, RenewableDispatch, RenewableFullDispatch)

    # define the renewable non-dispatch model
    PowerSimulations.set_device_model!(template_uc, RenewableNonDispatch, FixedOutput)
end

function define_device_model_weather_stochastic(template_uc; load_timeseries_name, renewable_timeseries_name)
    # Define custom PowerLoad model for MC runs
    load_model = DeviceModel(
        PowerLoad,
        StaticPowerLoad;
        time_series_names = Dict(ActivePowerTimeSeriesParameter => load_timeseries_name)
    )
    # Assign custom load model to template
    PowerSimulations.set_device_model!(template_uc, load_model)

    # Define custom RenewableDispatch model for MC runs (FTM)
    ftm_renewable_model = DeviceModel(
        RenewableDispatch,
        RenewableFullDispatch;
        time_series_names = Dict(ActivePowerTimeSeriesParameter => renewable_timeseries_name)
    )
    # Assign custom RenewableDispatch model to template
    PowerSimulations.set_device_model!(template_uc, ftm_renewable_model)

    # Define custom RenewableNonDispatch model for MC runs (BTM)
    btm_renewable_model = DeviceModel(
        RenewableNonDispatch,
        FixedOutput;
        time_series_names = Dict(ActivePowerTimeSeriesParameter => renewable_timeseries_name)
    )
    # Assign custom RenewableDispatch model to template
    PowerSimulations.set_device_model!(template_uc, btm_renewable_model)
end

function define_branch_model(template_uc)
    # define the branch model (static branch -> add unbounded flow variables and use flow constraints)
    PowerSimulations.set_device_model!(template_uc, AreaInterchange, StaticBranch)  #Q: is this defining TX?
end

function define_network_model(template_uc)
    # define our area balance power model to produce our zonal topology (Sienna default is copperplate)

    area_interchange = NetworkModel(
        AreaBalancePowerModel, #creates power balance constraints for each area
        use_slacks=false, #disable slack variables; i.e. prevent line flow exceedance
    )
    #assign the area_interchange object from above as our network model
    set_network_model!(template_uc, area_interchange)
end

function build_and_execute_simulation(template_uc, sys::System, paths::Dict; decision_name::String)
    UC_decision = DecisionModel(
        template_uc,
        sys;
        name = decision_name,
        optimizer = optimizer_with_attributes(Gurobi.Optimizer, "MIPGap" => 1e-2),
        system_to_file = false, # write the json and hf files
        initialize_model = true, # Q: what does this do?
        optimizer_solve_log_print = true, #solver output
        direct_mode_optimizer = true, # performance thing; default is true; set it false if you have specific need
        rebuild_model = false, # never have to use this, R&D thing
        store_variable_names = true,
        calculate_conflict = true, #infeasibility (gurobi only)
        export_optimization_model = false, # this exports the LP (location is...)
    )

    # Create Simulation Models
    sim_model = SimulationModels(
        decision_models = [UC_decision],
    )

    # Define Simulation Sequence
    sim_sequence = SimulationSequence(
        models = sim_model,
    )

    # Define Simulation
    sim = Simulation(
        name = "test-sim",
        steps = 3,  # Steps in your simulation
        models = sim_model,
        sequence = sim_sequence,
        simulation_folder = mktempdir(paths[:sim_files_dir], cleanup = true),
    )

    # Build the simulation folder
    build!(sim; console_level = Logging.Info,)

    # Execute the simulation
    execute!(sim, enable_progress_bar = true)

    # Return sim and UC_decision
    return sim, UC_decision

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

    # Export generation by fuel type
#=     all_gen_data = PowerAnalytics.get_generation_data(results) # do we need to define timeseries for all inputs for this to work?
    cat = PowerAnalytics.make_fuel_dictionary(sys)
    fuel = PowerAnalytics.categorize_data(all_gen_data.data, cat; curtailment = true, slacks = true)
    fuel_agg = PowerAnalytics.combine_categories(fuel)
    CSV.write(joinpath(path_scenario, "generation_by_fuel.csv"), fuel_agg) =#
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
    CSV.write(joinpath(paths[:data_dir], "nve_nameplate_capacity.csv"), df);
end
