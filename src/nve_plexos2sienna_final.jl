###########################
# Master Version of NVE-Sienna Script (single and multiple weather year)
###########################

# call packages
using Pkg
using Revise
using PowerSystems
using PowerSimulations
using PowerSystemCaseBuilder
using HydroPowerSimulations
using StorageSystemsSimulations
using InfrastructureSystems
using PowerAnalytics
using Gurobi # solver
using Dates
using Logging
using DataFrames
using Plots
using CSV
using TimeSeries
using SiennaPRASInterface

###########################
# Specify Run Type
###########################
run_type = "Monte_Carlo"; # "Deterministic" or "Monte_Carlo"

###########################
# Define Your Paths
###########################
# Specify the output name and scenario name dynamically
r2x_output_name = "output_stable_dec24"; # name of the r2x scenario
scenario_name = "output_stable_dec24";  # name of the sienna scenario

# Call the function to initialize paths and inputs
include(joinpath(@__DIR__, "NVE_non_weather.jl"));
paths = initialize_paths_and_inputs(r2x_output_name, scenario_name)

###########################
# Call the helpers function
###########################
# call _helpers.jl to load additional functions (modify_data, create_area_interchanges, attach_reserve_requirements, constrain_market_imports)
# define file path
helpers_file_path = joinpath(paths[:data_dir],"sienna_runs","_helpers.jl")
# load the script
include(helpers_file_path)

###########################
# Model Administration 
###########################
# Configure logging
logger = configure_logging(console_level=Logging.Info);

###########################
# Build The Initial System
###########################
sys, pw_data = build_system_from_r2x(paths, r2x_output_name);
get_units_base(sys)
set_units_base_system!(sys, "NATURAL_UNITS")

###########################
# Modify The Initial System
###########################
# create area interchanges, define initial SOC, set generator services, define availability, add VOM to production costs, set reference bus)
modify_system!(sys, pw_data, paths)

###########################
# Missing time series for fuel_prices for ThermalStandard objects 
###########################
# assign fuel_price timeseries to dual fuel Thermal standard objects (R2X defaulted to hydrogen fuel price)
active_unit = get_component(ThermalStandard, sys, "Valmy CT 3")
show_time_series(active_unit)
get_time_series_array(SingleTimeSeries, active_unit, "max_active_power")
get_time_series_array(SingleTimeSeries, active_unit, "fuel_price")
active_unit.operation_cost

# check to see if there are any missing timeseries for fuel_price
for g in collect(get_components(ThermalStandard, sys))
    if "fuel_price" ∉ get_name.(get_time_series_keys(g))
        @show g.name
        @show g.operation_cost.variable.fuel_cost
    end
end

# create and assign the timeseries for natural gas fuel prices
# ThermalStandard_missing_ts!(sys,paths)

# check to make sure missing timeseries were created successfully
#show_time_series(active_unit)
#get_time_series_array(SingleTimeSeries, active_unit, "fuel_price")
#active_unit.operation_cost

#= # remaining units should be oddballs (e.g., geothermal, waste heat, and biomass); assuming fixed fuel prices
for g in collect(get_components(ThermalStandard, sys))
    if "fuel_price" ∉ get_name.(get_time_series_keys(g))
        @show g.name
        #@show g.operation_cost.variable.fuel_cost
    end
end =#

###########################
# Convert DPV units from RenewableDispatch to RenewableNonDispatch
###########################
convert_BTM_units!(sys, paths)

###########################
# Define production costs for TA resource and add in DR Units
###########################
Tuning_Adjustment_Costs!(sys)
Demand_Response_CleanUp!(sys)

###########################
# Fix Hydro Dispatch Profile for Hoover
###########################
fix_Hydro_Dispatch!(sys)

###########################
# Spot Check your defined resources in the system
#= ###########################
# Initialize an empty DataFrame
df_gen = DataFrame(Resource = String[])
df_storage = DataFrame(Resource = String[])

# Loop through all generators in the system and collect their names
for generator in get_components(Generator, sys)
    name = get_name(generator)
    push!(df_gen, (name,))
end

for storage in get_components(Storage, sys)
    name = get_name(storage)
    push!(df_storage, (name,))
end

#combine the two dataframes
df_supply_stack = vcat(df_gen, df_storage)

# Write the DataFrame to a CSV file
CSV.write("NVE_supply_stack.csv", df_supply_stack) =#

############################################################
# Reassign prime mover type for ThermalStandard generators
############################################################
# read in datafile
df_pm = CSV.read(joinpath(paths[:data_dir],"nve_prime_mover_mapping.csv"), DataFrame);

# loop through all ThermalStandard generators
for thermal_gen in get_components(ThermalStandard, sys)
    obj_name = get_name(thermal_gen)  # Retrieve the name of the object

    # Find the corresponding row in the mapping table
    match_idx = findfirst(x -> x == obj_name, df_pm.resource)  # Match with the `resource` column
    if match_idx !== nothing  # case when a match is found
        # Access the prime_mover_type from the matched row
        new_prime_mover_type = df_pm.prime_mover_type[match_idx]  # Match with the `prime_mover_type` column
        set_prime_mover_type!(thermal_gen, new_prime_mover_type)
        # println("Updated $(obj_name) with PrimeMoverType: $new_prime_mover_type")  # Optional debug print
    else
        println("No match found for $(obj_name)")  # Debug print for unmatched objects
    end
end

#= #set all ThermalStandard generators that are CTs to zero
for thermal_gen in get_components(ThermalStandard, sys)
    if string(get_prime_mover_type(thermal_gen)) == "GT"
        set_active_power!(thermal_gen, 0.0)
    end
end =#

############################################################
# Assign Outage Statistics to all generators and storage devices
############################################################
# read in outage statistics
df_outage_stats = CSV.read(joinpath(paths[:LOLP_inputs],"outage_statistics.csv"), DataFrame)

#loop through generators
for generator in get_components(Generator, sys)
    # retrieve the object name 
    obj_name = get_name(generator)
    # Set the outage statistics
    outage_stats = df_outage_stats[df_outage_stats.Item .== obj_name, :]
    # check if outage_stats exist for active unit
    if size(outage_stats, 1) == 0
        @warn "No outage statistics found for $(obj_name)"
        continue
    end

    # define the outage transition data (i.e. MTTR and OTP)
    transition_data = GeometricDistributionForcedOutage(;
        mean_time_to_recovery=outage_stats.MTTR[1],  # Units of hours
        outage_transition_probability=outage_stats.outage_transition_probability[1],)  # Probability for outage per hour
    
    # Add the supplemental attribute to the generator
    add_supplemental_attribute!(sys, generator, transition_data)
end 

#loop through storage devices
for storage in get_components(Storage, sys)
    # retrieve the object name 
    obj_name = get_name(storage)
    # Set the outage statistics
    outage_stats = df_outage_stats[df_outage_stats.Item .== obj_name, :]
    # check if outage_stats is empty
    if size(outage_stats, 1) == 0
        @warn "No outage statistics found for $(obj_name)"
        continue
    end

    # define the outage transition data (i.e. MTTR and OTP)
    transition_data = GeometricDistributionForcedOutage(;
        mean_time_to_recovery=outage_stats.MTTR[1],  # Units of hours
        outage_transition_probability=outage_stats.outage_transition_probability[1],)  # Probability for outage per hour
    
    # Add the supplemental attribute to the generator
    add_supplemental_attribute!(sys, storage, transition_data)
end 

###########################
# Inspect Data
###########################
buses = collect(get_components(ACBus, sys));
loads = collect(get_components(PowerLoad, sys));
lines = collect(get_components(Line, sys));
arcs = collect(get_components(Arc, sys));
areas = collect(get_components(Area, sys));
area_interchanges = collect(get_components(AreaInterchange, sys));
gens_thermal = collect(get_components(ThermalStandard, sys));
gens_ftm_renew = collect(get_components(RenewableDispatch, sys));
gens_btm_renew = collect(get_components(RenewableNonDispatch, sys));
gen_hydro = collect(get_components(HydroDispatch, sys));
# solar_renew = [gen for gen in gens_renew if get_prime_mover_type(gen) == PrimeMovers.PVe];
# wind_renew = [gen for gen in gens_renew if get_prime_mover_type(gen) == PrimeMovers.WT];
# Initialize empty collections for solar and wind generators
solar_renew = []
wind_renew = []

# Loop through the renewable generators and classify them
for gen in gens_ftm_renew
    if get_prime_mover_type(gen) == PrimeMovers.PVe
        push!(solar_renew, gen)  # Add to solar_renew
    elseif get_prime_mover_type(gen) == PrimeMovers.WT
        push!(wind_renew, gen)  # Add to wind_renew
    end
end
batteries = collect(get_components(EnergyReservoirStorage, sys));
reserves_spinning = collect(get_components(VariableReserve{ReserveUp}, sys));
reserves_non_spinning = collect(get_components(VariableReserveNonSpinning, sys));

#geothermal inspection
active_unit = get_component(ThermalStandard, sys, "McGinness Hills Geothermal")
show_time_series(active_unit)
get_time_series_array(SingleTimeSeries, active_unit, "max_active_power")

# Retrieve all generator-type components in one go
all_generators = collect(get_components(Generator, sys))
all_storage = collect(get_components(Storage, sys))

# Define collections to iterate over
unit_collections = Dict(
    "GenUnits" => all_generators,
    "StorageUnits" => all_storage)

###########################
# Query Nameplate Capacity of System
###########################
system_capacity_query(unit_collections, paths);

############################################################
##  Timeseries
############################################################
show_time_series(sys)

#retrieve thermal standard objects with timeseries (fuel_price & max_active_power)
thermal_ts = collect(get_components(x -> has_time_series(x), ThermalStandard, sys));
active_unit = thermal_ts[1]
active_unit = get_component(ThermalStandard, sys, "Valmy CT 3")
show_time_series(active_unit)
get_time_series_keys(active_unit)
get_time_series_array(SingleTimeSeries, active_unit, "max_active_power"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_unit, "max_active_power"; ignore_scaling_factors = false)
get_time_series_array(SingleTimeSeries, active_unit, "fuel_price"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_unit, "fuel_price"; ignore_scaling_factors = false)

#retrieve RenewableDispatch objects with timeseries (rating factor, max_active_power)
renew_ts = collect(get_components(x -> has_time_series(x), RenewableDispatch, sys));
active_re_unit = get_component(RenewableDispatch, sys, "_PCM Generic Expansion_SPPC_PV")
active_re_unit = get_component(RenewableDispatch, sys, "Spring Valley Wind")
active_re_unit = get_component(RenewableDispatch, sys, "ACE Searchlight Solar")
show_time_series(active_re_unit)
get_time_series_keys(active_re_unit)
get_time_series_array(SingleTimeSeries, active_re_unit, "Rating Factor")
get_time_series_array(SingleTimeSeries, active_re_unit, "max_active_power"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_re_unit, "max_active_power"; ignore_scaling_factors = false)

#retrieve PowerLoad objects with timeseries
power_load_ts = collect(get_components(x -> has_time_series(x), PowerLoad, sys));
active_load = power_load_ts[1]
get_max_active_power(active_load)
show_time_series(active_load)
get_time_series_keys(active_load)
get_time_series_array(SingleTimeSeries, active_load, "max_active_power"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_load, "max_active_power"; ignore_scaling_factors = false)

active_load = power_load_ts[2]
get_max_active_power(active_load)
show_time_series(active_load)
get_time_series_keys(active_load)
get_time_series_array(SingleTimeSeries, active_load, "max_active_power"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_load, "max_active_power"; ignore_scaling_factors = false)

#= """
NOTE: base year  load fx - which comes from WY 1998 and is represented by "max_active_power" is slightly different
than what we will see with stochastic load fx, which is represented by "max_active_power_Y1998"

this reflects the same setup that was done in PLEXOS
""" =#

###########################
# Process TimeSeries / Forecasts 
###########################
if run_type == "Deterministic"
    # do nothing
elseif run_type == "Monte_Carlo"

    # define weather_related functions
    include("NVE_weather.jl")

    # FTM Solar
    ###########################
    #define solar PowerSystems.jl time series
    FTM_solar_ts_container, df_solar_max = define_FTM_solar_time_series(paths[:FTM_solar_dir], sys);

    # let's review our work
    df_FTM_solar = DataFrame(generator = map(x -> x[1], FTM_solar_ts_container),
                time_series = map(x -> x[2].name, FTM_solar_ts_container));

    # check to see if all objects are listed equal number of times
    df_FTM_solar_check = combine(groupby(df_FTM_solar, :generator), nrow => :count)

    #now we will add our solar time series container to the system
    add_FTM_solar_time_series_to_system!(sys,FTM_solar_ts_container);

    # now that we have the PowerSystems.jl timeseries added to the system, let's check your work
    active_re_unit = get_component(RenewableDispatch, sys, "_PCM Generic Expansion_SPPC_PV")
    show_time_series(active_re_unit)
    get_max_active_power(active_re_unit)
    get_time_series_array(SingleTimeSeries, active_re_unit, "Rating Factor"; ignore_scaling_factors = true)
    get_time_series_array(SingleTimeSeries, active_re_unit, "Rating Factor"; ignore_scaling_factors = false)
    get_time_series_array(SingleTimeSeries, active_re_unit, "max_active_power"; ignore_scaling_factors = true)
    get_time_series_array(SingleTimeSeries, active_re_unit, "max_active_power"; ignore_scaling_factors = false)
    get_time_series_array(SingleTimeSeries, active_re_unit, "max_active_power_Y1998"; ignore_scaling_factors = true)
    get_time_series_array(SingleTimeSeries, active_re_unit, "max_active_power_Y1998"; ignore_scaling_factors = false)

    # BTM Solar
    ###########################
    #define solar PowerSystems.jl time series
    BTM_solar_ts_container, df_BTM_solar_max = define_BTM_solar_time_series(paths[:BTM_solar_dir], sys);

    # let's review our work
    df_BTM_solar = DataFrame(generator = map(x -> x[1], BTM_solar_ts_container),
                time_series = map(x -> x[2].name, BTM_solar_ts_container));

    # check to see if all objects are listed equal number of times
    df_BTM_solar_check = combine(groupby(df_BTM_solar, :generator), nrow => :count)

    #now we will add our solar time series container to the system
    add_BTM_solar_time_series_to_system!(sys,BTM_solar_ts_container)

    # now that we have the PowerSystems.jl timeseries added to the system, let's check your work
    active_re_unit = get_component(RenewableNonDispatch, sys, "DPV_Nevada Power ND")
    show_time_series(active_re_unit)
    get_max_active_power(active_re_unit)
    get_time_series_array(SingleTimeSeries, active_re_unit, "Rating Factor"; ignore_scaling_factors = true)
    get_time_series_array(SingleTimeSeries, active_re_unit, "Rating Factor"; ignore_scaling_factors = false)
    get_time_series_array(SingleTimeSeries, active_re_unit, "max_active_power"; ignore_scaling_factors = true)
    get_time_series_array(SingleTimeSeries, active_re_unit, "max_active_power"; ignore_scaling_factors = false)
    get_time_series_array(SingleTimeSeries, active_re_unit, "max_active_power_Y1999"; ignore_scaling_factors = true)
    get_time_series_array(SingleTimeSeries, active_re_unit, "max_active_power_Y1999"; ignore_scaling_factors = false)

    # Wind
    ###########################
    wind_ts_container = define_wind_ts_container(paths[:wind_dir], sys);

    # let's review our work
    df_wind = DataFrame(generator = map(x -> x[1], wind_ts_container),
                time_series = map(x -> x[2].name, wind_ts_container));

    # check to see if all objects are listed equal number of times
    df_wind_check = combine(groupby(df_wind, :generator), nrow => :count)

    #now we will add our wind time series container to the system
    add_wind_time_series_to_system!(sys,wind_ts_container)

    #spot check your work
    # active_re_unit = get_component(RenewableDispatch, sys, "Spring Valley Wind")
    active_re_unit = get_component(RenewableDispatch, sys, "_PCM Generic Expansion_SPPC_Wind (ID)")
    show_time_series(active_re_unit)
    get_max_active_power(active_re_unit)
    #get_time_series_array(SingleTimeSeries, active_re_unit, "Rating Factor"; ignore_scaling_factors = true)
    #get_time_series_array(SingleTimeSeries, active_re_unit, "Rating Factor"; ignore_scaling_factors = false)
    get_time_series_array(SingleTimeSeries, active_re_unit, "max_active_power"; ignore_scaling_factors = true)
    get_time_series_array(SingleTimeSeries, active_re_unit, "max_active_power"; ignore_scaling_factors = false)
    get_time_series_array(SingleTimeSeries, active_re_unit, "max_active_power_Y1998"; ignore_scaling_factors = true)
    get_time_series_array(SingleTimeSeries, active_re_unit, "max_active_power_Y1998"; ignore_scaling_factors = false)

    #############
    # Demand
    #############
    #show PowerLoad objects w/ timeseries
    show_components(sys, PowerLoad, Dict("has_time_series" => x -> has_time_series(x)))

    #show time series for each PowerLoad object
    show_time_series(loads[1]) #default name "max_active_power"
    show_time_series(loads[2]) #default name "max_active_power"

    # create / upload powersystems.jl timeseries for regional areas of Nevada Power and Sierra regions
    load_ts_container = upload_split_load_forecasts(paths[:load_dir], loads);

    # check your work
    df_load = DataFrame(generator = map(x -> x[1], load_ts_container),
                time_series = map(x -> x[2].name, load_ts_container));

    # check to see if all objects are listed equal number of times
    df_PL_check = combine(groupby(df_load, :generator), nrow => :count)

    # Call the function to add time series to the system
    add_load_time_series_to_system!(sys, load_ts_container)

    #spot check your work
    active_load = get_component(PowerLoad, sys, "Nevada Power")
    show_time_series(active_load)
    get_max_active_power(active_load)
    get_time_series_array(SingleTimeSeries, active_load, "max_active_power"; ignore_scaling_factors = true)
    get_time_series_array(SingleTimeSeries, active_load, "max_active_power"; ignore_scaling_factors = false)
    get_time_series_array(SingleTimeSeries, active_load, "max_active_power_Y1998"; ignore_scaling_factors = true)
    get_time_series_array(SingleTimeSeries, active_load, "max_active_power_Y1998"; ignore_scaling_factors = false)

    active_load = get_component(PowerLoad, sys, "Sierra")
    show_time_series(active_load)
    get_max_active_power(active_load)
    get_time_series_array(SingleTimeSeries, active_load, "max_active_power"; ignore_scaling_factors = true)
    get_time_series_array(SingleTimeSeries, active_load, "max_active_power"; ignore_scaling_factors = false)
    get_time_series_array(SingleTimeSeries, active_load, "max_active_power_Y1998"; ignore_scaling_factors = true)
    get_time_series_array(SingleTimeSeries, active_load, "max_active_power_Y1998"; ignore_scaling_factors = false)

else
    @warn "Incorrect setting for run_type; $run_type is not a valid option"
end

#################################
# create timeseries fxs (48-hr horizon & 24-hr lookahead; i.e. 24 hour "realized intervals")
#################################
# create DeterministicSingleTimeSeries objects for all missing time series serving as forecasts 
transform_single_time_series!(sys, Hour(48), Hour(24)) 

###########################
# Sienna Temporary Patches
###########################
# first get baseline values
active_unit = get_component(ThermalStandard, sys, "Chuck Lenzi 2_B")
show_time_series(active_unit)
active_unit.operation_cost #note how fuel_cost has a fixed value specified (this is ignoring the ts we have attached)
show_time_series(active_load)

# correct for bug that has production cost expression ignoring fuel_prices sourced from timeseries
update_thermal_fuel_price_timeseries!(sys)

# check your work (ThermalStandard)
show_time_series(active_unit)
get_time_series(SingleTimeSeries, active_unit, "fuel_price")
get_time_series_array(SingleTimeSeries, active_unit, "fuel_price"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_unit, "fuel_price"; ignore_scaling_factors = false) # this should be equal to the previous cmd
get_time_series_array(DeterministicSingleTimeSeries, active_unit, "fuel_price")
active_unit.operation_cost # note: you should not see a fixed price under fuel_cost anymore; this should be a pointer to the timeseries

###########################################################################################
# For Loop to Iterate Over Each Weather Year for PCM in PowerSimualtions.jl WITHOUT FORCED OUTAGES
# This objective of running this portion of code is to ensure that all the weather-related parameters are flowing through correctly
###########################################################################################
# Define the range of weather years
if run_type == "Deterministic"
    weather_years = 1998;
elseif run_type == "Monte_Carlo" 
    weather_years = 1999:1999; # testing  only a few yrs to ensure proper configuration across weather years
else
    @warn "Incorrect setting for run_type; $run_type is not a valid option"
end 

# Iterate over each weather year
for wy in weather_years
    # Generate strings for the current year
    year_str = string(wy)

    if run_type == "Deterministic"

        #assign name
        uc_decision_name = "base_yr"

        # Create an empty model reference
        template_uc = ProblemTemplate()

        # Define thermal, hydro, and storage device models (non-weather-dependent)
        define_device_model_non_weather(template_uc)

        # Define load and renewable device models (deterministic)
        define_device_model_weather_deterministic(template_uc)

        # Define branch and network models
        define_branch_model(template_uc)
        define_network_model(template_uc)
    elseif run_type == "Monte_Carlo"
        #define unique identifies for timeseries pointers and decision model names
        uc_decision_name = "DA_WY_$year_str"  # e.g., "DA_WY_1998"
        ts_WY_name = "max_active_power_Y$year_str"  # e.g., "max_active_power_Y1998"

        #################################
        # Define Device, Branch, an Network Models
        #################################
        # Create an empty model reference
        template_uc = ProblemTemplate()

        # Define thermal, hydro, and storage device models (non-weather-dependent)
        define_device_model_non_weather(template_uc)

        # Define load and renewable dispatch models - both FTm and BTM (weather-dependent)
        define_device_model_weather_stochastic(template_uc; load_timeseries_name=ts_WY_name, renewable_timeseries_name=ts_WY_name)

        # Define branch and network models
        define_branch_model(template_uc)
        define_network_model(template_uc)
    else
        @warn "Incorrect setting for run_type; $run_type is not a valid option"
    end

    ###########################
    # Build and Execute Simulation
    ###########################
    get_units_base(sys)
    set_units_base_system!(sys, "NATURAL_UNITS")
    sim, UC_decision = build_and_execute_simulation(template_uc, sys, paths; decision_name=uc_decision_name);

    # Print a message to indicate that the simulation is complete
    println("Simulation completed for: $uc_decision_name")

    ###########################
    # Export the Results
    ###########################
    if run_type == "Deterministic"
        # define file paths to store (processed) results for the active year
        file_path = (paths[:scenario_dir_d])
    elseif run_type == "Monte_Carlo"
        # define file paths to store (processed) results for the active weather yr
        file_path = joinpath(paths[:scenario_dir_s], "results_WY_$year_str")
    else
        @warn "Incorrect setting for run_type; $run_type is not a valid option"
    end 
    
    # check if folder directory exists; if not, create it
    if !ispath(file_path)
        mkpath(file_path)
    else
        # do nothing
    end

    # export the results
    query_write_export_results(sim, file_path, uc_decision_name)        
end

###########################################################################################
# For Loop to Iterate Over Each Weather Year for our Monte Carlo RA LOLP Simulation using SiennaPRASInterface
###########################################################################################
# Define the range of weather years
weather_years = 1998:2022 # testing for only a few yrs for now
# Dictionary to store results for each weather year
# weather_year_results = Dict{Int, Tuple{Simulation, DecisionModel}}()
shortfall_results = []
shortfall_stats_array = []

# Iterate over each weather year
for wy in weather_years

    # manual entry (i.e. for troubleshooting purposes)
    #wy = 1998

    # Generate strings for the current year
    year_str = string(wy)

    #define unique identifies for timeseries pointers and decision model names
    uc_decision_name = "DA_WY_$year_str"  # e.g., "DA_WY_1998"
    ts_WY_name = "max_active_power_Y$year_str"  # e.g., "max_active_power_Yxxxx"

    ###########################
    # Run SiennaPRASInterface
    ###########################
    # DeviceRAModel (PowerSystems Device, SPI Abstract Formulation, time_series_names, kwargs)
    problem_template = RATemplate(
        PowerSystems.Area,
        [
            #define thermal device model
            DeviceRAModel(PowerSystems.ThermalStandard,GeneratorPRAS(max_active_power="max_active_power"),),

            #define hydro device model
            DeviceRAModel(PowerSystems.HydroGen, GeneratorPRAS(max_active_power="max_active_power"),),

            #define storage device model
            DeviceRAModel(PowerSystems.EnergyReservoirStorage, EnergyReservoirLossless()),
                
            #define renewable device model (both RenewableDispatch and RenewableNonDispatch) - weather dependent
            DeviceRAModel(PowerSystems.RenewableGen, GeneratorPRAS(max_active_power=ts_WY_name),),

            #define load device model (weather dependent)
            DeviceRAModel(PowerSystems.PowerLoad, StaticLoadPRAS(max_active_power=ts_WY_name),),

            #define line device model
            DeviceRAModel(PowerSystems.Line, LinePRAS(),),

            #define AreaInterchange device model
            DeviceRAModel(PowerSystems.AreaInterchange, AreaInterchangeLimit(),),
        ],
    )    

    # Convert PowerSimulations.jl system (i.e., PSY) to PRAS model (i.e., PSI) 
    pras_sys = generate_pras_system(sys, problem_template)

    ##########################
    #SPI Troubleshooting (make sure load is being sent over to PRAS in NATURAL UNITS)
    ##########################
    const SPI = SiennaPRASInterface
    const PSY = PowerSystems
    loads_to_formula = SPI.build_component_to_formulation(SPI.LoadPRAS, sys, problem_template.device_models)
    static_ts_summary = PSY.get_static_time_series_summary_table(sys)
    s2p_meta = SPI.S2P_metadata(static_ts_summary)
    regions= get_components(Area, sys)
    SPI.get_region_loads(s2p_meta,regions, loads_to_formula,)

    pras_sys.generators.names
    # find index (idx) of generator of interest in pras_sys.generators
    idx = findfirst(==("DPV_Nevada Power ND"), pras_sys.generators.names)
    idx = findfirst(==("_PCM Generic Expansion_SPPC_Wind (ID)"), pras_sys.generators.names)
    profile = pras_sys.generators.capacity[idx,:]

    active_re_unit = get_component(RenewableNonDispatch, sys, "DPV_Nevada Power ND")    
    active_re_unit = get_component(RenewableDispatch, sys, "_PCM Generic Expansion_SPPC_Wind (ID)")
    get_max_active_power(active_re_unit)
    psy_no_scaling = values(get_time_series_array(SingleTimeSeries, active_re_unit, "max_active_power_Y1999"; ignore_scaling_factors = true))
    psy_scaling = values(get_time_series_array(SingleTimeSeries, active_re_unit, "max_active_power_Y1999"; ignore_scaling_factors = false))

    df_summary = DataFrame(pras_ts = profile, psy_no_scaling = psy_no_scaling, psy_scaling = psy_scaling) 

    CSV.write("generator_capacity.csv", df_summary)

    #General Inquiry (appliestodevice checks)
    typeof(pras_sys)
    DeviceRAModel_type = typeof(DeviceRAModel(PSY.ThermalStandard,GeneratorPRAS(),))   

    SPI.appliestodevice(DeviceRAModel(PSY.ThermalStandard,GeneratorPRAS(),), typeof(active_re_unit))
    SPI.appliestodevice(DeviceRAModel(PSY.ThermalStandard,GeneratorPRAS(),), typeof(active_unit))

    SPI.appliestodevice(DeviceRAModel(PSY.Generator,GeneratorPRAS(),), typeof(active_re_unit))
    SPI.appliestodevice(DeviceRAModel(PSY.Generator,GeneratorPRAS(),), typeof(active_unit))

    #define the PRAS method -> sequential MC with 20 samples of 8760 FO draws
    # Set Sequential Monte Carlo method as method along with # of samples and seed value
    # note: Sequential MC honors chronology
    method = SequentialMonteCarlo(samples=20, seed=1, verbose=true, threaded=false)

    # Run the PRAS simulation and assess shortfall periods and related statistics
    shortfalls, shortfall_stats = assess(pras_sys, method, ShortfallSamples(), Shortfall())
    push!(shortfall_results, shortfalls)
    push!(shortfall_stats_array, shortfall_stats)
    println("SiennaPRASInterface simulation completed for: WY $wy.")
end

###################################
#Export Results for post-processing 
###################################
# Create a DataFrame to store summary statistics for all weather years
summary_df = DataFrame(weather_year = Int[], eue = Float64[], std_error = Float64[])
hourly_ens = DataFrame(hour = 1:8760)

# Process each weather year's results
for (idx, wy) in enumerate(weather_years)
    # Calculate EUE (Expected Unserved Energy) for this weather year
    eue_value = EUE(shortfall_stats_array[idx])
    eue = eue_value.eue.estimate
    std_error = eue_value.eue.standarderror
    push!(summary_df, (wy, eue, std_error))

    # Sum along the first dimension
    summed_shortfall = dropdims(sum(shortfall_results[idx].shortfall, dims=1), dims=1) # 8760 by nsamples array

    # Create a new DataFrame for this weather year's data
    year_df = DataFrame(summed_shortfall, :auto)

    # Rename columns with weather year prefix
    colnames = [Symbol("$(wy)_x$i") for i in 1:size(summed_shortfall, 2)]
    DataFrames.rename!(year_df, colnames)

    # Horizontally concatenate with the main DataFrame
    hourly_ens = hcat(hourly_ens, year_df)
end

#####################
# write Shortfall Summary across all weather years to csv
#####################
summary_csv_path = joinpath(paths[:scenario_dir_s], "shortfall_summary.csv")
CSV.write(summary_csv_path, summary_df)
println("Exported shortfall summary statistics to: $(summary_csv_path)")

#####################
# write hourly ENS across all samples to csv
#####################
# Write the combined DataFrame to CSV
all_samples_csv_path = joinpath(paths[:scenario_dir_s], "all_shortfall_samples.csv")
CSV.write(all_samples_csv_path, hourly_ens)

#####################
# calculate daily MWh and export to csv
#####################
# Drop the :hour column if it's just 1 to 8760 and not needed
df_daily_ENS = select(hourly_ens, Not(:hour))
# Add a new column for day number (1 to 365)
df_daily_ENS.Day = ceil.(Int, (1:nrow(df_daily_ENS)) ./ 24)
# Group by Day and sum each group
df_daily_ENS = combine(groupby(df_daily_ENS, :Day), names(df_daily_ENS, Not(:Day)) .=> sum)
# Create a Date column for the year 2030
df_daily_ENS.Date = Date(2030, 1, 1) .+ Day.(df_daily_ENS.Day .- 1)
#reorder columns
select!(df_daily_ENS, :Date, Not([:Date]))
#drop Day
select!(df_daily_ENS, Not([:Day]))
#export the daily shortfall amt
daily_shortfall_csv_path = joinpath(paths[:scenario_dir_s], "daily_shortfall_ENS.csv")
CSV.write(daily_shortfall_csv_path, df_daily_ENS)

#####################
# calculate daily ENS hours and export to csv
#####################
# Drop the :hour column if it's just 1 to 8760 and not needed
df_daily_hrs = select(hourly_ens, Not(:hour))
# Add a new column for day number (1 to 365)
df_daily_hrs.Day = ceil.(Int, (1:nrow(df_daily_hrs)) ./ 24)
# Group by Day and sum each group
df_daily_hrs = combine(groupby(df_daily_hrs, :Day), names(df_daily_hrs, Not(:Day)) .=> x -> count(>(0), x))
# Create a Date column for the year 2030
df_daily_hrs.Date = Date(2030, 1, 1) .+ Day.(df_daily_hrs.Day .- 1)
#reorder columns
select!(df_daily_hrs, :Date, Not([:Date]))
#drop Day
select!(df_daily_hrs, Not([:Day]))
#export the daily shortfall amt
daily_shortfall_csv_path = joinpath(paths[:scenario_dir_s], "daily_shortfall_hours.csv")
CSV.write(daily_shortfall_csv_path, df_daily_hrs)

#####################
# calculate daily event indicators and export to csv
#####################
# Drop the :hour column if it's just 1 to 8760 and not needed
df_daily_events = select(hourly_ens, Not(:hour))
# Add a new column for day number (1 to 365)
df_daily_events.Day = ceil.(Int, (1:nrow(df_daily_events)) ./ 24)
df_daily_events
# Group by Day and sum each group
df_daily_events = combine(groupby(df_daily_events, :Day),names(df_daily_events, Not(:Day)) .=> (x -> Int(any(>(0), x))))
# Create a Date column for the year 2030
df_daily_events.Date = Date(2030, 1, 1) .+ Day.(df_daily_events.Day .- 1)
#reorder columns
select!(df_daily_events, :Date, Not([:Date]))
#drop Day
select!(df_daily_events, Not([:Day]))
#export the daily shortfall amt
daily_shortfall_csv_path = joinpath(paths[:scenario_dir_s], "daily_shortfall_hours_indicator.csv")
CSV.write(daily_shortfall_csv_path, df_daily_events)


