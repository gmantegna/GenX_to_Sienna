module GenX_to_Sienna_FPA

using PowerSystems
using PowerSimulations
using PowerAnalytics
using Gurobi
using CSV
using DataFrames
using TimeSeries
using Dates
using Statistics
using StorageSystemsSimulations
using HydroPowerSimulations
using DataStructures 

##########################
# Include helper functions
##########################
include(joinpath(@__DIR__, "helpers.jl"))
include(joinpath(@__DIR__, "system_builder.jl"))
include(joinpath(@__DIR__, "time_series_builder.jl"))
include(joinpath(@__DIR__, "simulation_builder.jl"))

##########################
# Define Paths
##########################
paths = initialize_paths_and_inputs()

##########################
# Initialize System
##########################
# initialize system with base power of 100MVA (aides in per unit calculations)
sys = System(100.0)
sys_base_power = get_base_power(sys)
get_units_base(sys)
set_units_base_system!(sys, "NATURAL_UNITS")
get_units_base(sys)

##########################
# Define Network Topology 
##########################
# GENX Network Topology
network_df = CSV.read(joinpath(paths[:data_dir], "system", "Network.csv"), DataFrame)

# Buses
##########################
network_mapping_df = network_df[:, [first(names(network_df)), "zone_num"]]
#rename the first column to "zone"
DataFrames.rename!(network_mapping_df, :Column1 => :zone)
network_mapping_df = dropmissing(network_mapping_df)
zone_dict = OrderedDict{String,Int}(
    string(row.zone) => row.zone_num 
    for row in eachrow(network_mapping_df))

# define PSY bus objects
buses_dict = create_buses(zone_dict)

# add buses to system
for (zone_name, bus) in buses_dict
    add_component!(sys, bus)
end

# define collection of buses
buses = collect(get_components(ACBus, sys))

# let's check our work
get_components(ACBus, sys)
show_components(ACBus, sys)
active_component = get_component(ACBus, sys, "BANC")
get_number(active_component)

# Areas
##########################
# let's first create our dictionary of area objects
areas_dict = create_areas(buses_dict)

# now we can add the areas to system and assign buses to areas
for (area_name, area) in areas_dict
    # add area to system
    add_component!(sys, area)
    
    # assign bus to area
    bus = get_component(ACBus, sys, area_name)
    set_area!(bus, area)
end

# let's check our work
get_components(Area, sys)
show_components(Area, sys)
active_component = collect(get_components(Area, sys))[1]
get_name(active_component)

# let's verify bus to area assignments
active_component = collect(get_components(ACBus, sys))[2]
get_name(get_area(active_component))

# define collection of areas
areas = collect(get_components(Area, sys))

# Lines
##########################
# pull in lines info from GENX Network Topology
existing_lines_mapping_df = network_df[:, Between("Network_Lines", "Profile_Reverse")]

# Check to see if Network expansion is active (i.e. new incremental transfer service is available)
if isfile(joinpath(paths[:data_dir], "results", "network_expansion.csv"))
    network_expansion=true
else
    network_expansion=false
end

if network_expansion
    candidate_lines_mapping_df = CSV.read(joinpath(paths[:data_dir], "results", "network_expansion.csv"), DataFrame);
else
    candidate_lines_mapping_df = nothing
end

# define PSY line  objects
lines_dict = create_lines(existing_lines_mapping_df, candidate_lines_mapping_df, sys)

# add lines to system (note: arc count should be equal to line count)
for (line_number, line) in lines_dict
    add_component!(sys, line)
end

# let's check our work
get_components(Line, sys)
show_components(Line, sys)
active_component = collect(get_components(Line, sys))[1]
get_name(active_component)
get_arc(active_component)
get_rating(active_component)

# define collection of lines
lines = collect(get_components(Line, sys));

# Area Interchanges
##########################
# create area interchanges
area_interchanges_dict = create_area_interchanges(existing_lines_mapping_df, sys_base_power, lines_dict)

# add area interchanges to system
for (area_interchange_name, area_interchange) in area_interchanges_dict
    add_component!(sys, area_interchange)
end

#= # remove all area interchanges from system
for Area_Interchange in collect(get_components(AreaInterchange, sys))
    remove_component!(sys, Area_Interchange)
end =#

# let's check our work
get_components(AreaInterchange, sys)
show_components(AreaInterchange, sys)
active_component = collect(get_components(AreaInterchange, sys))[1]
get_name(active_component)
get_from_area(active_component)
get_to_area(active_component)
get_flow_limits(active_component)

# define collection of area interchanges
area_interchanges = collect(get_components(AreaInterchange, sys))

# Transmission Interfaces
##########################
# read interface data
interface_df = CSV.read(joinpath(paths[:data_dir], "system", "Simultaneous_Flow_Constraints.csv"), DataFrame)

# create transmission interfaces dictionary
interfaces_dict = create_transmission_interfaces(interface_df, lines_dict, sys_base_power)

# add interfaces to system
for (interface_name, interface) in interfaces_dict
    add_component!(sys, interface)
end

# let's check our work
get_components(TransmissionInterface, sys)
show_components(TransmissionInterface, sys)
active_component = collect(get_components(TransmissionInterface, sys))[1]
get_name(active_component)
get_active_power_flow_limits(active_component)
get_direction_mapping(active_component)

# define collection of transmission interfaces
transmission_interfaces = collect(get_components(TransmissionInterface, sys))

##########################
# Define PowerLoads 
##########################
# read in relevant CSV files

# define file_path 
demand_data_path = joinpath(paths[:data_dir], "system", "Demand_data.csv")
# call the function to generate the demand timeseries df
demand_ts = process_demand_data(demand_data_path, zone_dict)

# Let's create our power loads dictionary
power_loads_dict = create_power_loads(demand_ts, sys)

# Now let's add our PowerLoad objects to the system
for (pl_name, pl_object) in power_loads_dict
    add_component!(sys, pl_object)
end

# Let's check our work
get_components(PowerLoad, sys)
show_components(PowerLoad, sys)
active_component = get_component(PowerLoad, sys, "Load_PGE")
get_base_power(active_component)

##########################
# Define Fuel Objects 
##########################
# define df of fuel mapping
fuel_mapping_df = CSV.read(joinpath(paths[:data_dir], "FuelMapping.csv"), DataFrame)

##########################
# Define Generators 
##########################

# General Gen-Related Info 
##########################
# retrieve capacity data
capacity_df = CSV.read(joinpath(paths[:data_dir], "results", "capacity.csv"), DataFrame);
# read in prime mover types mapping
PM_type_df = CSV.read(joinpath(paths[:data_dir], "MoverTypesMapping.csv"), DataFrame);
# create prime mover type dictionary
PM_type_dict = Dict((row.Key) => row.Value for row in eachrow(PM_type_df));

# create storage type dictionary
storage_type_df = CSV.read(joinpath(paths[:data_dir], "StorageMapping.csv"), DataFrame);
storage_type_dict = Dict((row.Key) => row.Value for row in eachrow(storage_type_df));

# Thermal Generators 
##########################
# read in thermal data
thermal_df = CSV.read(joinpath(paths[:data_dir], "resources", "Thermal.csv"), DataFrame);

# define thermal generator objects
ThermalStandard_dict = create_ThermalStandard_objects(sys, thermal_df, capacity_df, PM_type_dict, fuel_mapping_df, zone_dict)
    
# add thermal generators to system
for (thermal_name, thermal_object) in ThermalStandard_dict
    add_component!(sys, thermal_object)
end

# define collection of thermal generators
ThermalStandard_generators = collect(get_components(ThermalStandard, sys))

# Check your work
active_component = ThermalStandard_generators[5]

get_name(active_component)
get_base_power(active_component) #installed nameplate capacity (MW)   
show_time_series(active_component) # no time series attached to this component (yet)
active_component.operation_cost #note how fuel_cost has a fixed value specified (this is ignoring the ts we have attached)
show_time_series(active_component)
# get_time_series(DeterministicSingleTimeSeries, active_component, "fuel_price")
# get_time_series_array(DeterministicSingleTimeSeries, active_component, "fuel_price")


# Renewable Dispatch Generators (i.e., VRE) 
##########################
# read in renewable  data
vre_df = CSV.read(joinpath(paths[:data_dir], "resources", "Vre.csv"), DataFrame);

# define Vre generator objects
Vre_dict = create_VRE_objects(sys, vre_df, capacity_df, PM_type_dict, zone_dict)

# add VRE generators to system
for (vre_name, vre_object) in Vre_dict
    add_component!(sys, vre_object)
end

# define collection of VRE generators
Renew_D_generators = collect(get_components(RenewableDispatch, sys))

# Check your work
active_component = Renew_D_generators[5]

get_name(active_component)
get_base_power(active_component) #installed nameplate capacity (MW)   
show_time_series(active_component) # no time series attached to this component (yet)
active_component.operation_cost #note how fuel_cost has a fixed value specified (this is ignoring the ts we have attached)
show_time_series(active_component) # no time series attached to this component (yet)

# Renewable NonDispatch Generators (i.e., BTM) 
##########################
# Check to see if Must run is active (must run is BTM PV)
if isfile(joinpath(paths[:data_dir], "resources", "Must_run.csv"))

    btm_df = CSV.read(joinpath(paths[:data_dir], "resources", "Must_run.csv"), DataFrame);

    # define Vre generator objects
    btm_dict = create_btm_objects(sys, btm_df, capacity_df, PM_type_dict, zone_dict)

    # add VRE generators to system
    for (btm_name, btm_object) in btm_dict
        add_component!(sys, btm_object)
    end

    # define collection of VRE generators
    Renew_ND_generators = collect(get_components(RenewableNonDispatch, sys))

    # Check your work
    active_component = Renew_ND_generators[3]

    get_name(active_component)
    get_base_power(active_component) #installed nameplate capacity (MW)   
    show_time_series(active_component) # no time series attached to this component (yet)

else
    # do nothing
end

# Hydro Generators 
##########################
# read in hydro generator data
hydro_df = CSV.read(joinpath(paths[:data_dir], "resources", "Hydro.csv"), DataFrame);

# define Vre generator objects
Hydro_Dispatch_dict = create_Hydro_objects(sys, hydro_df, capacity_df, PM_type_dict, zone_dict)

# add VRE generators to system
for (hydro_name, hydro_object) in Hydro_Dispatch_dict
    add_component!(sys, hydro_object)
end

# define collection of Hydro  generators
HydroDispatch_generators = collect(get_components(HydroDispatch, sys))

# Check your work
active_component = HydroDispatch_generators[3]

get_name(active_component)
get_base_power(active_component) #installed nameplate capacity (MW)   
show_time_series(active_component) # no time series attached to this component (yet)
active_component.operation_cost #note how fuel_cost has a fixed value specified (this is ignoring the ts we have attached)
show_time_series(active_component) # no time series attached to this component (yet)

# Storage Resources 
##########################
# read in storage data      
storage_df = CSV.read(joinpath(paths[:data_dir], "resources", "Storage.csv"), DataFrame);

# define storage generator objects
Storage_dict = create_storage_objects(sys, storage_df, capacity_df, PM_type_dict, storage_type_dict, zone_dict)    

# add storage generators to system
for (storage_name, storage_object) in Storage_dict
    add_component!(sys, storage_object)
end

# define collection of storage generators
Storage_objects = collect(get_components(Storage, sys));

# Check your work
active_component = Storage_objects[5]    

get_name(active_component)
get_base_power(active_component) #installed nameplate capacity (MW)   
show_time_series(active_component) # no time series attached to this component (yet)
active_component.operation_cost #note how fuel_cost has a fixed value specified (this is ignoring the ts we have attached)


##########################
# Define time series for PSY objects
##########################

# General
##########################
# define file_path 
generator_variability_data_path = joinpath(paths[:data_dir], "system", "Generators_variability.csv")

# call the function to generate the demand timeseries df
gen_variability_df = process_generator_variability_data(generator_variability_data_path)

# let's write the gen_variability_df to a csv file
CSV.write(joinpath(paths[:data_dir], "Generators_variability.csv"), gen_variability_df)

# Power Loads
##########################
# first let's create our PSI timeseries objects and store them in a container structured as a nested dictionary
PL_ts_container = create_demand_PSY_timeseries(demand_ts, power_loads_dict)

# Now we add those PSY timeseries to the PowerLoad objects in the system
for (device_name, year_ts_dict) in PL_ts_container
    # Retrieve the active device by its name
    active_device = get_component(PowerLoad, sys, "Load_"*device_name)

    if active_device !== nothing
        # Loop through each year's time series for this device
        for (year, time_series) in year_ts_dict
            # Add the time series to the system
            add_time_series!(sys, active_device, time_series)
            println("Added time series: ", time_series.name, " for year ", year, " to device: Load_", device_name)
        end
    else
        @warn "Device $device_name not found in the system. Time series not added."
    end
end

# let's check our work
active_load = get_component(PowerLoad, sys, "Load_PGE")
show_time_series(active_load) # now we have time series attached to the PowerLoad object
get_time_series_array(SingleTimeSeries, active_load, "max_active_power_1998"; ignore_scaling_factors = true) #p.u.; units: device base 
get_time_series_array(SingleTimeSeries, active_load, "max_active_power_1998"; ignore_scaling_factors = false) # natural units

# Renewable Dispatch Generators
##########################
# Create dictionary of time series for renewable dispatch generators
renewable_ts_container = create_Renew_D_PSY_timeseries(gen_variability_df, Renew_D_generators)

# spot check the time series
active_ts = renewable_ts_container["Idaho_Wind_PGE"]["1998"]
active_ts.name
active_ts.data


# Add all the time series to system
for (resource_name, year_ts_dict) in renewable_ts_container
    # Retrieve the active device by its name
    active_device = get_component(RenewableDispatch, sys, resource_name)

    if active_device !== nothing
        # Loop through each year's time series for this resource
        for (year, time_series) in year_ts_dict
            # Add the time series to the system
            add_time_series!(sys, active_device, time_series)
            println("Added time series: ", time_series.name, " for year ", year, " to device: ", resource_name)
        end
    else
        @warn "Device $resource_name not found in the system. Time series not added."
    end
end

# Let's check our work
active_object = Renew_D_generators[1]
show_time_series(active_object)
ts_key = get_time_series_keys(active_object)
ts_ref = get_time_series_keys(active_object).ref
ts_size = get_time_series_keys(active_object).size
get_time_series_array(SingleTimeSeries, active_object, "max_active_power_1998"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_object, "max_active_power_1998"; ignore_scaling_factors = false)

#= # remove time series of RenewableDispatch objects from system
for gen in Renew_D_objects # loop through the collection of RenewableDispatch objects
    for i in length(get_time_series_keys(gen))
        # retrieve the time series
        ts_key = get_time_series_keys(gen)[i]
        ts_name = get_name(ts_key)
        # remove the time series
        remove_time_series!(sys, SingleTimeSeries, gen, ts_name)
    end
end =#

# Renewable Non-Dispatch Generators
##########################
# Create dictionary of time series for renewable non-dispatch generators
renew_ND_ts_container = create_Renew_ND_PSY_timeseries(gen_variability_df,Renew_ND_generators)

# spot check the time series
active_ts = renew_ND_ts_container["Customer_PV_PGE"]["1998"]
active_ts.name
active_ts.data

# Add all the time series to system
for (resource_name, year_ts_dict) in renew_ND_ts_container
    # Retrieve the active device by its name
    active_device = get_component(RenewableNonDispatch, sys, resource_name)

    if active_device !== nothing
        # Loop through each year's time series for this resource
        for (year, time_series) in year_ts_dict
            # Add the time series to the system
            add_time_series!(sys, active_device, time_series)
            println("Added time series: ", time_series.name, " for year ", year, " to device: ", resource_name)
        end
    else
        @warn "Device $resource_name not found in the system. Time series not added."
    end
end

# Let's check our work
active_object = Renew_ND_generators[1]
show_time_series(active_object)
ts_key = get_time_series_keys(active_object)
ts_ref = get_time_series_keys(active_object).ref
ts_size = get_time_series_keys(active_object).size
get_time_series_array(SingleTimeSeries, active_object, "max_active_power_1998"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_object, "max_active_power_1998"; ignore_scaling_factors = false)

# Thermal Standard Generators
##########################
# Create dictionary of time series for thermal standard non-dispatch generators
ThermalStandard_ts_container = create_ThermalStandard_PSY_timeseries(gen_variability_df,ThermalStandard_generators)

# spot check the time series
active_ts = ThermalStandard_ts_container["CAISO_Aero_CT_PGE"]["1998"]
active_ts.name
active_ts.data

# Add all the time series to system
for (resource_name, year_ts_dict) in ThermalStandard_ts_container
    # Retrieve the active device by its name
    active_device = get_component(ThermalStandard, sys, resource_name)

    if active_device !== nothing
        # Loop through each year's time series for this resource
        for (year, time_series) in year_ts_dict
            # Add the time series to the system
            add_time_series!(sys, active_device, time_series)
            println("Added time series: ", time_series.name, " for year ", year, " to device: ", resource_name)
        end
    else
        @warn "Device $resource_name not found in the system. Time series not added."
    end
end

# Let's check our work
active_object = ThermalStandard_generators[1]
show_time_series(active_object)
ts_key = get_time_series_keys(active_object)
ts_ref = get_time_series_keys(active_object).ref
ts_size = get_time_series_keys(active_object).size
get_time_series_array(SingleTimeSeries, active_object, "max_active_power_1998"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_object, "max_active_power_1998"; ignore_scaling_factors = false)

# Hydro Generators
##########################
# Create dictionary of time series for thermal standard non-dispatch generators
Hydro_ts_container = create_Hydro_PSY_timeseries(gen_variability_df,HydroDispatch_generators)

# spot check the time series
active_ts = Hydro_ts_container["CAISO_Hydro_PGE"]["1998"]
active_ts.name
active_ts.data

# Add all the time series to system
for (resource_name, year_ts_dict) in Hydro_ts_container
    # Retrieve the active device by its name
    active_device = get_component(HydroDispatch, sys, resource_name)

    if active_device !== nothing
        # Loop through each year's time series for this resource
        for (year, time_series) in year_ts_dict
            # Add the time series to the system
            add_time_series!(sys, active_device, time_series)
            println("Added time series: ", time_series.name, " for year ", year, " to device: ", resource_name)
        end
    else
        @warn "Device $resource_name not found in the system. Time series not added."
    end
end

# Let's check our work
active_object = HydroDispatch_generators[1]
show_time_series(active_object)
ts_key = get_time_series_keys(active_object)
ts_ref = get_time_series_keys(active_object).ref
ts_size = get_time_series_keys(active_object).size
get_time_series_array(SingleTimeSeries, active_object, "max_active_power_1998"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_object, "max_active_power_1998"; ignore_scaling_factors = false)

# remove time series (nuclear option)
#remove_time_series!(sys, SingleTimeSeries)


#################################
# create timeseries fxs 
#################################
# create DeterministicSingleTimeSeries objects (48-hr horizon & 24-hr lookahead; i.e. 24 hour "realized intervals") 
transform_single_time_series!(sys, Hour(48), Hour(24)) 

##########################
# Define Fuel Info 
##########################
# define df of fuel prices 
fuels_df = CSV.read(joinpath(paths[:data_dir], "system", "Fuels_data.csv"), DataFrame)

# call the function to generate the demand timeseries df
# fuel_ts = process_fuel_data(fuels_df)


##########################
# Define Outage Data
##########################

##########################
# Define Reserves
##########################
# general
reserves_df = CSV.read(joinpath(paths[:data_dir], "system", "Operational_reserves.csv"), DataFrame);

# Spinning Reserves
spin_requirement = reserves_df[!, "Rsv_Req_Percent_Demand"];

# Reserve Zones
oprsv_zones = CSV.read(joinpath(paths[:data_dir], "oprsv_zones.csv"), DataFrame).Zone; #CAISO TAC Zones

##########################
# Define PowerSimulations.jl (PSI) template and model 
##########################
# define run_type
run_type = "Deterministic"

# determine if run_type is deterministic or monte-create
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
    wy = 1998

    if run_type == "Deterministic"

        #assign name
        uc_decision_name = "deterministic_$wy"

        # Create an empty model reference
        template_uc = ProblemTemplate()

        # Define non-weather related Device Models
        ##########################
        # storage
        define_storage_model(template_uc)

        # Define weather-dependent Device Models
        ##########################
        # thermal
        define_thermal_model(template_uc, wy)

        # hydro
        define_hydro_model(template_uc, wy)

        # load
        define_load_model(template_uc, wy)

        # renewable dispatch
        define_renewable_dispatch_model(template_uc, wy)

        # renewable non-dispatch
        define_renewable_non_dispatch_model(template_uc, wy)

        # Define branch model
        define_branch_model(template_uc)

        # Define network model
        define_network_model(template_uc)

    elseif run_type == "Monte_Carlo"
        # do nothing
        # TO-DO: layer in Sienna-PRAS Interface
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

end # module
