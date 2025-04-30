#################################
# Define Device Models in PSI 
#################################
function define_storage_model(template_uc)
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

function define_PHS_model(template_uc)
    # Define custom DeviceModel for EnergyReservoirStorage
    PHS_model = DeviceModel(
        HydroPumpedStorage,
        StorageDispatchWithReserves,
    )
    # Assign the storage model to the template_uc
    PowerSimulations.set_device_model!(template_uc, PHS_model)
end


function define_thermal_model(template_uc, WY)
    # Define custom DeviceModel for ThermalStandard with time-varying max_active_power
#=     thermal_model = DeviceModel(ThermalStandard, ThermalStandardUnitCommitment;
        time_series_names = Dict(PowerSimulations.ActivePowerTimeSeriesParameter => "max_active_power_$WY")) =#

    thermal_model = DeviceModel(ThermalStandard, ThermalStandardUnitCommitment; time_series_names = Dict{Any, String}(
                PowerSimulations.FuelCostParameter => "fuel_price",
                PowerSimulations.ActivePowerTimeSeriesParameter => "max_active_power_$WY",))

    # assign the thermal model to the template_uc
    PowerSimulations.set_device_model!(template_uc, thermal_model)
end

function define_hydro_model(template_uc, WY)
    # Define Hydro model
    hydro_model = DeviceModel(HydroDispatch,HydroDispatchRunOfRiver;
        time_series_names = Dict(PowerSimulations.ActivePowerTimeSeriesParameter => "max_active_power_$WY"))
    
    # assign the hydro model to the template_uc
    PowerSimulations.set_device_model!(template_uc, hydro_model)
end

function define_load_model(template_uc, WY)
    # define the load model with time series
    load_model = DeviceModel(PowerLoad, StaticPowerLoad;
        time_series_names = Dict(PowerSimulations.ActivePowerTimeSeriesParameter => "max_active_power_$WY"))
    # Assign load model to template
    PowerSimulations.set_device_model!(template_uc, load_model)
end

function define_renewable_dispatch_model(template_uc, WY)
    # define the renewable dispatch model with time series
    renewable_dispatch_model = DeviceModel(RenewableDispatch, RenewableFullDispatch;
        time_series_names = Dict(PowerSimulations.ActivePowerTimeSeriesParameter => "max_active_power_$WY"))

    # Assign renewable dispatch model to template
    PowerSimulations.set_device_model!(template_uc, renewable_dispatch_model)
end

function define_renewable_non_dispatch_model(template_uc, WY)
    # define the renewable non-dispatch model with time series
    renewable_non_dispatch_model = DeviceModel(RenewableNonDispatch, FixedOutput;
        time_series_names = Dict(PowerSimulations.ActivePowerTimeSeriesParameter => "max_active_power_$WY"))

    # Assign renewable non-dispatch model to template
    PowerSimulations.set_device_model!(template_uc, renewable_non_dispatch_model)
end

#################################
# Define Branch Model in PSI 
#################################
function define_branch_model(template_uc)
    # define the branch model to be assigned to our AreaInterchanges
    # static branch -> adds unbounded flow variables and uses flow constraints

    # Q: how would i add slack variables to this?
    AI_branch = DeviceModel(AreaInterchange, StaticBranch, use_slacks=false)

    PowerSimulations.set_device_model!(template_uc, AI_branch)
end



#################################
# Define Network Model in PSI 
#################################
function define_AreaNetwork_model(template_uc)
    # define our area balance power model to honor our zonal topology bc default in Sienna is copperplate

    Zonal_AI = NetworkModel(
        AreaBalancePowerModel, #Approximation to represent inter-area flow with each area represented as a single node.
        use_slacks=false, #disable slack variables; i.e. prevent line flow exceedance
    )
    #assign the area_interchange object from above as our network model
    set_network_model!(template_uc, Zonal_AI)
end

function define_CopperPlate_model(template_uc)
    # define our area balance power model to honor our zonal topology bc default in Sienna is copperplate

    CopperSheet = NetworkModel(
        CopperPlatePowerModel, #Approximation to represent inter-area flow with each area represented as a single node.
        use_slacks=true, 
    )
    #assign the area_interchange object from above as our network model
    set_network_model!(template_uc, CopperSheet)
end



#################################
# Define Service Models in PSI 
#################################
# regulation up
function define_RegUp_service_model(template_uc)
    # Define the regulation up reserve service model with time series requirements
    reg_reserve_up_model = ServiceModel(
        PSY.VariableReserve{PSY.ReserveUp},
        RangeReserve;
        time_series_names = Dict(PowerSimulations.RequirementTimeSeriesParameter => "requirement")
    )

    # Assign the service model to the template_uc
    PowerSimulations.set_service_model!(template_uc, reg_reserve_up_model)
end

# regulation down
function define_RegDown_service_model(template_uc)
    # Define the regulation down reserve service model with time series requirements
    reg_reserve_down_model = ServiceModel(
        VariableReserve{ReserveDown},
        RangeReserve;
        time_series_names = Dict(PowerSimulations.RequirementTimeSeriesParameter => "requirement")
    )

    # Assign the service model to the template_uc
    PowerSimulations.set_service_model!(template_uc, reg_reserve_down_model)
end


#################################
# Define Simulation Model in PSI 
#################################
function build_and_execute_simulation(template_uc::ProblemTemplate, sys::System, paths::Dict, decision_name::String)
    # initialize our decision model
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

    # Initialize Simulation Model(s)
    sim_model = SimulationModels(
        decision_models = [UC_decision],
    )

    # Initialize Simulation Sequence
    sim_sequence = SimulationSequence(
        models = sim_model,
    )

    # Define the simulation
    sim = Simulation(
        name = "test-sim",
        steps = 3,  # Steps in your simulation
        models = sim_model,
        sequence = sim_sequence,
        simulation_folder = mktempdir(paths[:output_dir_base], cleanup = true),
    )

    # Build the simulation folder
    build!(sim; console_level = Logging.Info,)

    # Execute the simulation
    execute!(sim, enable_progress_bar = true)

    # Return sim and UC_decision
    return sim, UC_decision

end
