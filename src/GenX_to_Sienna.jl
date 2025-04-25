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
using Logging

##########################
# Include helper functions
##########################
include(joinpath(@__DIR__, "helpers.jl"))
include(joinpath(@__DIR__, "system_builder.jl"))
include(joinpath(@__DIR__, "time_series_builder.jl"))
include(joinpath(@__DIR__, "simulation_builder.jl"))

##########################
# Model Administration
##########################
# Define Path Directories
paths = initialize_paths_and_inputs()

# Define logger
logger = configure_logging(console_level=Logging.Info);

##########################
# Initialize System
##########################
# initialize system with base power of 100MVA (aides in per unit calculations)
sys = System(1.0) # 100 is the default but we are using 1.0 bc defining in natural units
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
# define collection of lines
lines = collect(get_components(Line, sys));

active_component = lines[1]
get_name(active_component)
get_arc(active_component)
get_rating(active_component)

# Create and write line parameters to CSV
file_path = joinpath(paths[:data_dir], "line_parameters.csv")
create_line_parameters_df(lines, file_path)

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
# define collection of area interchanges
area_interchanges = collect(get_components(AreaInterchange, sys))

active_component = area_interchanges[1]
get_name(active_component)
get_from_area(active_component)
get_to_area(active_component)
get_flow_limits(active_component)
get_flow_limits(active_component).from_to

# Create and write area interchange parameters to CSV
file_path = joinpath(paths[:data_dir], "area_interchange_parameters.csv")
create_area_interchange_parameters_df(area_interchanges, file_path)

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
# define collection of transmission interfaces
transmission_interfaces = collect(get_components(TransmissionInterface, sys))
active_component = transmission_interfaces[1]
get_name(active_component)
get_active_power_flow_limits(active_component)
get_direction_mapping(active_component)

# Create and write transmission interface parameters to CSV
file_path = joinpath(paths[:data_dir], "transmission_interface_parameters.csv");
create_transmission_interface_parameters_df(sys, paths);

##########################
# Define PowerLoads 
##########################
# read in relevant CSV files

# define file_path 
demand_data_path = joinpath(paths[:data_dir], "system", "Demand_data.csv")
# call the function to generate the demand timeseries df
demand_ts = process_demand_data(demand_data_path, zone_dict)

# Let's create our power loads dictionary
power_loads_dict = create_power_loads(demand_ts, sys);

# Now let's add our PowerLoad objects to the system
for (pl_name, pl_object) in power_loads_dict
    add_component!(sys, pl_object)
end

# Let's check our work
get_components(PowerLoad, sys)
show_components(PowerLoad, sys)
active_component = get_component(PowerLoad, sys, "Load_PGE")
get_name(active_component)
get_base_power(active_component)
get_active_power(active_component)
get_bus(active_component)

# define collection of power loads
power_loads = collect(get_components(PowerLoad, sys))

# Create and write power load parameters to CSV
file_path = joinpath(paths[:data_dir], "powerload_parameters.csv");
create_powerload_parameters_df(sys, paths);


# remove all powerloads from system
#= for Power_Load in collect(get_components(PowerLoad, sys))
    remove_component!(sys, Power_Load)
end =#

#= # let's write the power loads to a csv file
CSV.write(joinpath(paths[:data_dir], "Power_loads.csv"), demand_ts) =#

##########################
# Define Fuel Objects 
##########################
# define df of fuel mapping
fuel_mapping_df = CSV.read(joinpath(paths[:data_dir], "FuelMapping.csv"), DataFrame)

##########################
# Define Generators / Storage Devices 
##########################

# General Info 
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
ThermalStandard_dict = create_ThermalStandard_objects(sys, thermal_df, capacity_df, PM_type_dict, fuel_mapping_df, zone_dict);
    
# add thermal generators to system
for (thermal_name, thermal_object) in ThermalStandard_dict
    add_component!(sys, thermal_object)
end

#= # remove all thermal generators from system
for Thermal_Standard in collect(get_components(ThermalStandard, sys))
    remove_component!(sys, Thermal_Standard)
end =#

# define collection of thermal generators
ThermalStandard_generators = collect(get_components(ThermalStandard, sys))

# Check your work
active_component = ThermalStandard_generators[5]

get_name(active_component)
get_base_power(active_component) #installed nameplate capacity (MW)   
get_bus(active_component)
get_rating(active_component)
get_active_power_limits(active_component)
get_time_limits(active_component)
get_fuel(active_component)
active_component.operation_cost #note how fuel_cost has a fixed value specified (this is ignoring the ts we have attached)

show_time_series(active_component)
# get_time_series(DeterministicSingleTimeSeries, active_component, "fuel_price")
# get_time_series_array(DeterministicSingleTimeSeries, active_component, "fuel_price")

#check to make sure no units with base power of 0 are in the system
show_components(ThermalStandard, sys, [:base_power])

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

#check to make sure no units with base power of 0 are in the system
show_components(RenewableDispatch, sys, [:base_power])

#= # remove all RenewableDispatch from system
for Renewable_Dispatch in collect(get_components(RenewableDispatch, sys))
    remove_component!(sys, Renewable_Dispatch)
end =#

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

#check to make sure no units with base power of 0 are in the system
show_components(RenewableNonDispatch, sys, [:base_power])

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

# check to make sure no units with base power of 0 are in the system
show_components(HydroDispatch, sys, [:base_power])

# remove all HydroDispatch from system
#= for Hydro_Dispatch in collect(get_components(HydroDispatch, sys))
    remove_component!(sys, Hydro_Dispatch)
end =#

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
Storage_objects = collect(get_components(EnergyReservoirStorage, sys));

# Check your work
active_component = Storage_objects[5]
get_name(active_component)
get_base_power(active_component) #installed nameplate capacity (MW) 
get_storage_capacity(active_component)    
get_storage_level_limits(active_component)
get_initial_storage_capacity_level(active_component)
get_efficiency(active_component)
get_input_active_power_limits(active_component)
get_output_active_power_limits(active_component)
get_rating(active_component)
active_component.operation_cost # check operation cost

# check to make sure no units with base power of 0 are in the system
show_components(EnergyReservoirStorage, sys, [:base_power])

# remove all storage objects from system
#= for Storage in collect(get_components(Storage, sys))
    remove_component!(sys, Storage)
end =#

# Troubleshooting: Create DataFrame with storage parameters and write to CSV
file_path = joinpath(paths[:data_dir], "storage_parameters.csv");
create_storage_parameters_df(Storage_objects, file_path);

# Pumped Hydro Storage Resources 
##########################

# define storage generator objects
PumpedHydro_dict = create_pumped_hydro_objects(sys, storage_df, capacity_df, PM_type_dict, zone_dict)   

# add storage generators to system
for (PHS_name, PHS_object) in PumpedHydro_dict
    add_component!(sys, PHS_object)
end

# define collection of storage generators
PHS_objects = collect(get_components(HydroPumpedStorage, sys));

# Check your work
active_component = PHS_objects[5]
get_name(active_component)
get_base_power(active_component) #installed nameplate capacity (MW) 
get_storage_capacity(active_component)
get_initial_storage(active_component)
get_pump_efficiency(active_component)
get_active_power_limits(active_component)
get_rating(active_component)
active_component.operation_cost # check operation cost

# check to make sure no units with base power of 0 are in the system
show_components(HydroPumpedStorage, sys, [:base_power])

# remove all PHS_objects from system
#= for PHS in collect(get_components(HydroPumpedStorage, sys))
    remove_component!(sys, PHS)
end =#

# Troubleshooting: Create DataFrame with storage parameters and write to CSV
file_path = joinpath(paths[:data_dir], "PHS_parameters.csv");
create_pumped_hydro_parameters_df(PHS_objects, file_path);


###########################
# Query Nameplate Capacity of System
###########################
# Retrieve all generator-type components in one go
all_generators = collect(get_components(Generator, sys));
all_storage = collect(get_components(Storage, sys));

# Define collections to iterate over
unit_collections = Dict(
    "GenUnits" => all_generators,
    "StorageUnits" => all_storage)

system_capacity_query(unit_collections, paths);

##########################
# Define time series for PSY objects
##########################

# General
##########################
# define file_path 
generator_variability_data_path = joinpath(paths[:data_dir], "system", "Generators_variability.csv");

# call the function to generate the generator profile timeseries df
gen_variability_df = process_generator_variability_data(generator_variability_data_path)

#= # let's write the gen_variability_df to a csv file
CSV.write(joinpath(paths[:data_dir], "Generators_variability.csv"), gen_variability_df) =#

# Power Loads
##########################
# first let's create our PSI timeseries objects and store them in a container structured as a nested dictionary
PL_ts_container = create_demand_PSY_timeseries(demand_ts, power_loads_dict)

#spot check the time series
active_ts = PL_ts_container["PGE"]["1998"]
active_ts.name
active_ts.data

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
active_ts = renewable_ts_container["Southern_NV_Eldorado_Solar_SCE"]["1998"]
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
active_object = Renew_D_generators[5]
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
active_ts = ThermalStandard_ts_container["CAISO_CCGT1_PGE"]["1998"]
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
active_object = ThermalStandard_generators[5]
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

# Pumped Hydro  Generators
##########################
# GENX has this pinned at 1 for all timesteps across all weather years; therefore not assigning any time series to the PHS objects


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

wy = weather_years

#assign name
decision_name = "deterministic_$wy"

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
###########################
define_branch_model(template_uc)

# Define network model
###########################
# CopperPlate
# define_CopperPlate_model(template_uc)

# AreaInterchange
define_AreaNetwork_model(template_uc)


#################################
# Define Simulation Model in PSI 
#################################
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
    simulation_folder = mktempdir(paths[:sienna_simulation_dir], cleanup = true),
)

# Build the simulation folder
build!(sim; console_level = Logging.Info,) # this will give us a "built" build status; run status still "initialized"

# Execute the simulation
execute!(sim, enable_progress_bar = true) # run status will now be "successfully_finalized" if successful

# Print a message to indicate that the simulation is complete
println("Simulation completed for: $decision_name")

###########################
# Export the Results
###########################
results_file_path = joinpath(paths[:sienna_results_dir], "results_$wy")
simulation_file_path = paths[:sienna_simulation_dir]


# check if results folder directory exists; if not, create it
if !ispath(results_file_path)
    mkpath(results_file_path)
else
    # do nothing
end

# check if simulation folder directory exists; if not, create it
if !ispath(simulation_file_path)
    mkpath(simulation_file_path)
else
    # do nothing
end

# export the results
# query_write_export_results(sim, file_path, uc_decision_name)

###########################
# Query Results
###########################
sim_results = SimulationResults(sim)
results = get_decision_problem_results(sim_results, decision_name) 

# Input TimeSeries Parameters
load_parameter = read_realized_parameter(results, "ActivePowerTimeSeriesParameter__PowerLoad")
thermal_parameter = read_realized_parameter(results, "ActivePowerTimeSeriesParameter__ThermalStandard")
renewDispatch_parameter = read_realized_parameter(results, "ActivePowerTimeSeriesParameter__RenewableDispatch")
renewNonDispatch_parameter = read_realized_parameter(results, "ActivePowerTimeSeriesParameter__RenewableNonDispatch") # Sienna doesnt store nonDispatch
hydro_parameter = read_realized_parameter(results, "ActivePowerTimeSeriesParameter__HydroDispatch")

# Output Realized Generation Values
thermal_active_power = read_realized_variable(results, "ActivePowerVariable__ThermalStandard")
renewDispatch_active_power = read_realized_variable(results, "ActivePowerVariable__RenewableDispatch")
hydro_active_power = read_realized_variable(results, "ActivePowerVariable__HydroDispatch")
storage_charge = read_realized_variable(results, "ActivePowerInVariable__EnergyReservoirStorage")
storage_discharge = read_realized_variable(results, "ActivePowerOutVariable__EnergyReservoirStorage")


# combine all FTM generators
gen_active_power = hcat(thermal_active_power, select(renewDispatch_active_power, Not(1)), select(hydro_active_power, Not(1)))

# Output Realized TX flows
AreaInterchange_flow = read_realized_variable(results, "FlowActivePowerVariable__AreaInterchange")
# TO-DO: Activate interfaces for SFCs

# Output Expressions
power_balance = read_realized_expression(results, "ActivePowerBalance__Area")

# Get Production Costs
pc_thermal = read_realized_expression(results, "ProductionCostExpression__ThermalStandard")
pc_renewable = read_realized_expression(results, "ProductionCostExpression__RenewableDispatch")
pc_hydro = read_realized_expression(results, "ProductionCostExpression__HydroDispatch")
pc_all = hcat(pc_thermal,select(pc_renewable, Not(1)), select(pc_hydro, Not(1)))

###########################
# Export Results
###########################
# Define output paths and write dataframes to CSV
CSV.write(joinpath(results_file_path, "load_active_power.csv"), load_parameter); # Input time series values
CSV.write(joinpath(results_file_path, "thermal_parameters.csv"), thermal_parameter); # Input time series values
CSV.write(joinpath(results_file_path, "FTM_renewable_parameters.csv"), renewDispatch_parameter); # Input time series values
CSV.write(joinpath(results_file_path, "BTM_active_power.csv"), renewNonDispatch_parameter); # Input time series values
CSV.write(joinpath(results_file_path, "hydro_parameter.csv"), hydro_parameter); # Input time series values

CSV.write(joinpath(results_file_path, "FTM_generator_active_power.csv"), gen_active_power);
CSV.write(joinpath(results_file_path, "storage_charge.csv"), storage_charge);
CSV.write(joinpath(results_file_path, "storage_discharge.csv"), storage_discharge);
CSV.write(joinpath(results_file_path, "AreaInterchange_flow.csv"), AreaInterchange_flow);
CSV.write(joinpath(results_file_path, "power_balance.csv"), power_balance);
CSV.write(joinpath(results_file_path, "production_costs.csv"), all_pc);    

    

end # module
