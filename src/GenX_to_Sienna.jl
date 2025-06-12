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
using InfrastructureSystems

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
# define constants 
const PSY = PowerSystems
const PSI = PowerSimulations
#const SPI = SiennaPRASInterface

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
bus_region_ba_df = CSV.read(joinpath(paths[:data_dir], "BusRegionBAMapping.csv"), DataFrame)

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
active_object = get_component(ACBus, sys, "BANC")
get_number(active_object)

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
active_object = collect(get_components(Area, sys))[1]
get_name(active_object)

# let's verify bus to area assignments
active_object = collect(get_components(ACBus, sys))[2]
get_name(get_area(active_object))

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

active_object = lines[1]
get_name(active_object)
get_arc(active_object)
get_rating(active_object)

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

active_object = area_interchanges[1]
get_name(active_object)
get_from_area(active_object)
get_to_area(active_object)
get_flow_limits(active_object)
get_flow_limits(active_object).from_to

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
active_object = transmission_interfaces[1]
get_name(active_object)
get_active_power_flow_limits(active_object)
get_direction_mapping(active_object)

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
demand_ts_df = process_demand_data(demand_data_path, zone_dict)

# let's write the power loads to a csv file
CSV.write(joinpath(paths[:data_dir], "demand_ts_df_8760.csv"), demand_ts_df)

# Let's create our power loads dictionary
power_loads_dict = create_power_loads(demand_ts_df, sys);

# Now let's add our PowerLoad objects to the system
for (pl_name, pl_object) in power_loads_dict
    add_component!(sys, pl_object)
end

# Let's check our work
get_components(PowerLoad, sys)
show_components(PowerLoad, sys)
active_object = get_component(PowerLoad, sys, "Load_PGE")
get_name(active_object)
get_base_power(active_object)
get_active_power(active_object)
get_bus(active_object)

# define collection of power loads
power_loads = collect(get_components(PowerLoad, sys));

# Create and write power load parameters to CSV
file_path = joinpath(paths[:data_dir], "powerload_parameters.csv");
create_powerload_parameters_df(sys, paths);

# remove all powerloads from system
#= for Power_Load in collect(get_components(PowerLoad, sys))
    remove_component!(sys, Power_Load)
end =#

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
# Read fuel mapping information and convert to dictionary
fuel_mapping_df = CSV.read(joinpath(paths[:data_dir], "FuelMapping.csv"), DataFrame)

# create storage type dictionary
storage_type_df = CSV.read(joinpath(paths[:data_dir], "StorageMapping.csv"), DataFrame);
storage_type_dict = Dict((row.Key) => row.Value for row in eachrow(storage_type_df));

# Thermal Generators 
##########################
# read in thermal data
thermal_df = CSV.read(joinpath(paths[:data_dir], "resources", "Thermal.csv"), DataFrame);

# Update fuel type for geothermal resources
for i in 1:nrow(thermal_df)
    if occursin("Geothermal", thermal_df[i, :Resource])
        thermal_df[i, :Fuel] = "Geothermal"
    end
end

# define thermal generator objects
ThermalStandard_dict = create_ThermalStandard_objects(sys, thermal_df, capacity_df, PM_type_dict, fuel_mapping_df, zone_dict);
    
# add thermal generators to system
for (thermal_name, thermal_object) in ThermalStandard_dict
    add_component!(sys, thermal_object)
end

# define collection of thermal generators
ThermalStandard_generators = collect(get_components(ThermalStandard, sys));

active_object = get_component(ThermalStandard, sys, "CAISO_Biomass_PGE")

# Update must-run status for baseload ThermalStandard objects (biomass, geothermal, nuclear, CHP)
update_TS_must_run_status!(ThermalStandard_generators)

# Check your work by visually inspecting the must_run status of the thermal generators
show_components(ThermalStandard, sys, [:must_run])

# check other attributes of a ThermalStandard
active_object = get_component(ThermalStandard, sys, "CAISO_CCGT1_PGE")

get_name(active_object)
get_base_power(active_object) #installed nameplate capacity (MW)   
get_bus(active_object)
get_rating(active_object)
get_active_power_limits(active_object)
get_time_limits(active_object)
get_fuel(active_object)
active_object.operation_cost #note how fuel_cost has a fixed value specified (this is ignoring the ts we have attached)

show_time_series(active_object)
# get_time_series(DeterministicSingleTimeSeries, active_object, "fuel_price")
# get_time_series_array(DeterministicSingleTimeSeries, active_object, "fuel_price")

#check to make sure no units with base power of 0 are in the system
show_components(ThermalStandard, sys, [:base_power])

# remove all thermal generators from system
#= for Thermal_Standard in collect(get_components(ThermalStandard, sys))
    remove_component!(sys, Thermal_Standard)
end =#


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
active_object = Renew_D_generators[5]

get_name(active_object)
get_base_power(active_object) #installed nameplate capacity (MW)   
show_time_series(active_object) # no time series attached to this component (yet)
active_object.operation_cost #note how fuel_cost has a fixed value specified (this is ignoring the ts we have attached)
show_time_series(active_object) # no time series attached to this component (yet)

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
    active_object = Renew_ND_generators[3]

    get_name(active_object)
    get_base_power(active_object) #installed nameplate capacity (MW)   
    show_time_series(active_object) # no time series attached to this component (yet)

else
    # do nothing
end

#check to make sure no units with base power of 0 are in the system
show_components(RenewableNonDispatch, sys, [:base_power])

# Hydro Generators (HydroDispatch) 
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
active_object = HydroDispatch_generators[1]

get_name(active_object)
get_base_power(active_object) #installed nameplate capacity (MW)   
show_time_series(active_object) # no time series attached to this component (yet)
active_object.operation_cost #note how fuel_cost has a fixed value specified (this is ignoring the ts we have attached)
show_time_series(active_object) # no time series attached to this component (yet)

#= # Check for missing hydro_budget time series and handle existing ones
for hydro in HydroDispatch_generators
    println(get_name(hydro))
    # show_time_series(hydro)
    if !has_time_series(hydro)
        @warn "Hydro resource $(get_name(hydro)) has no timeseries"
    elseif has_time_series(hydro)
        #check if "hydro_budget" timeseries is in the system
        if "hydro_budget" ∈ get_name.(get_time_series_keys(hydro))
            println("hydro_budget time series found")
        else
            println("hydro_budget time series not found")
        end
    end
end =#

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
active_object = Storage_objects[5]
get_name(active_object)
get_base_power(active_object) #installed nameplate capacity (MW) 
get_storage_capacity(active_object)    
get_storage_level_limits(active_object)
get_initial_storage_capacity_level(active_object)
get_efficiency(active_object)
get_input_active_power_limits(active_object)
get_output_active_power_limits(active_object)
get_rating(active_object)
active_object.operation_cost # check operation cost

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
PumpedHydro_dict = create_PHS_objects(sys, storage_df, capacity_df, PM_type_dict, zone_dict)   

# add storage generators to system
for (PHS_name, PHS_object) in PumpedHydro_dict
    add_component!(sys, PHS_object)
end

# define collection of storage generators
PHS_objects = collect(get_components(HydroPumpedStorage, sys));

# Check your work
active_object = PHS_objects[5]
get_name(active_object)
get_base_power(active_object) #installed nameplate capacity (MW) 
get_storage_capacity(active_object)
get_initial_storage(active_object)
get_pump_efficiency(active_object)
get_active_power_limits(active_object)
get_rating(active_object)
get_status(active_object)
# set_status!(active_object, PSY.PumpHydroStatusModule.PumpHydroStatus.GEN)
#set_status!(active_object, PSY.PumpHydroStatusModule.PumpHydroStatus.OFF)
#get_status(active_object)
active_object.operation_cost # check operation cost
# check to make sure no units with base power of 0 are in the system
show_components(HydroPumpedStorage, sys, [:base_power])

# remove all PHS_objects from system
#= for PHS in collect(get_components(HydroPumpedStorage, sys))
    remove_component!(sys, PHS)
end =#

# Troubleshooting: Create DataFrame with storage parameters and write to CSV
file_path = joinpath(paths[:data_dir], "PHS_parameters.csv");
create_PHS_parameters_df(PHS_objects, file_path);

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
gen_variability_ts_df = process_generator_variability_data(generator_variability_data_path)


#= # TROUBLESHOOTING:  the gen_variability_ts_df to a csv file
CSV.write(joinpath(paths[:data_dir], "Generators_variability.csv"), gen_variability_ts_df) =#

# Power Loads
##########################
# first let's create our PSI timeseries objects and store them in a container structured as a nested dictionary
PL_ts_container = create_demand_PSY_timeseries(demand_ts_df, power_loads_dict)

#spot check the time series
active_ts = PL_ts_container["PGE"]["2001"]
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
get_time_series_array(SingleTimeSeries, active_load, "max_active_power_2001"; ignore_scaling_factors = true) #p.u.; units: device base 
get_time_series_array(SingleTimeSeries, active_load, "max_active_power_2001"; ignore_scaling_factors = false) # natural units

# Renewable Dispatch Generators
##########################
# Create dictionary of time series for renewable dispatch generators
renewable_ts_container = create_Renew_D_PSY_timeseries(gen_variability_ts_df, Renew_D_generators)

# spot check the time series
active_ts = renewable_ts_container["Southern_NV_Eldorado_Solar_SCE"]["2001"]
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
active_object = get_component(RenewableDispatch, sys, "CAISO_Solar_PGE")
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
renew_ND_ts_container = create_Renew_ND_PSY_timeseries(gen_variability_ts_df,Renew_ND_generators)

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
get_time_series_array(SingleTimeSeries, active_object, "max_active_power_2001"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_object, "max_active_power_2001"; ignore_scaling_factors = false)

# Thermal Standard Generators
##########################
# Create dictionary of time series for thermal standard non-dispatch generators
ThermalStandard_ts_container = create_ThermalStandard_PSY_timeseries(gen_variability_ts_df,ThermalStandard_generators)

# spot check the time series
active_ts = ThermalStandard_ts_container["CAISO_CCGT1_PGE"]["2001"]
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
get_time_series_array(SingleTimeSeries, active_object, "max_active_power_2001"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_object, "max_active_power_2001"; ignore_scaling_factors = false)

# Hydro Max Active Power
##########################

# Create dictionary of time series for thermal standard non-dispatch generators
Hydro_ts_container = create_Hydro_PSY_timeseries(gen_variability_ts_df,HydroDispatch_generators)

# spot check the time series
active_ts = Hydro_ts_container["CAISO_Hydro_PGE"]["2001"]
active_ts.name
active_ts.data

# Add all the max_active_power time series to system
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

# Hydro Energy Budgets
##########################
# define budget df timeseries
hydro_budget_ts_df = process_hydro_budget_data(joinpath(paths[:data_dir], "system", "Hourly_energy_budget.csv"))

# Create dictionary of time series for thermal standard non-dispatch generators
Hydro_budget_ts_container = create_Hydro_Budget_PSY_timeseries(hydro_budget_ts_df,HydroDispatch_generators)

# spot check the time series
active_ts = Hydro_budget_ts_container["CAISO_Hydro_PGE"]["2001"]
active_ts.name
active_ts.data

# Add all the hydro budget time series to system
for (resource_name, year_ts_dict) in Hydro_budget_ts_container 
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
active_object = get_component(HydroDispatch, sys, "CAISO_Hydro_PGE")
show_time_series(active_object)
ts_key = get_time_series_keys(active_object)
ts_ref = get_time_series_keys(active_object).ref
ts_size = get_time_series_keys(active_object).size
get_time_series_array(SingleTimeSeries, active_object, "max_active_power_1998"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_object, "max_active_power_1998"; ignore_scaling_factors = false)

# Pumped Hydro  Generators
##########################
# GENX has max_active_power  pinned at 1 for all timesteps across all weather years;
# however we do need to create our PSY timeseries for the inflows and outflows

# Create dictionary of time series for pumped hydro generators
PHS_ts_container = create_PHS_PSY_timeseries(PHS_objects)

# spot check the time series
active_ts = PHS_ts_container["CAISO_Pumped_Hydro_PGE"]["inflow"]
active_ts.name
active_ts.data

active_ts = PHS_ts_container["CAISO_Pumped_Hydro_PGE"]["outflow"]
active_ts.name
active_ts.data

# Add all inflow and outflow time series to PHS objects
for (resource_name, year_ts_dict) in PHS_ts_container
    # Retrieve the active device by its name
    active_device = get_component(HydroPumpedStorage, sys, resource_name)

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

# check our work
active_object = get_component(HydroPumpedStorage, sys, "CAISO_Pumped_Hydro_PGE")
show_time_series(active_object)
ts_key = get_time_series_keys(active_object)
get_time_series_array(SingleTimeSeries, active_object, "inflow"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_object, "outflow"; ignore_scaling_factors = true)

#= # remove time series from PHS objects in system
for gen in PHS_objects # loop through the collection of PHS objects
    for i in length(get_time_series_keys(gen))
        # retrieve the time series
        ts_key = get_time_series_keys(gen)[i]
        ts_name = get_name(ts_key)
        # remove the time series
        remove_time_series!(sys, SingleTimeSeries, gen, ts_name)
    end
end =#

# now let's define our inflow and outflow timeseries for the PHS objects

##########################
# Define Reserves
##########################
# general
##########################
reserves_df = CSV.read(joinpath(paths[:data_dir], "system", "Operational_reserves.csv"), DataFrame);
oprsv_zones_df = CSV.read(joinpath(paths[:data_dir], "oprsv_zones.csv"), DataFrame)
oprsv_zones = oprsv_zones_df.Zone

# Regulation
##########################
# intialize services  
reg_reserve_serv_dict = create_CAISO_reg_reserve_services()
reg_reserve_serv_dict["CAISO_reg_up"]
reg_reserve_serv_dict["CAISO_reg_down"]

# define contributing devices
reg_reserve_units_dict = create_CAISO_reg_reserve_units(ThermalStandard_generators, Renew_D_generators, HydroDispatch_generators, Storage_objects, PHS_objects, thermal_df, vre_df, hydro_df, storage_df)
reg_reserve_units_dict["CAISO_reg_up"]
reg_reserve_units_dict["CAISO_reg_down"]

# add reg services to system
add_service!(sys, reg_reserve_serv_dict["CAISO_reg_up"], reg_reserve_units_dict["CAISO_reg_up"])
add_service!(sys, reg_reserve_serv_dict["CAISO_reg_down"], reg_reserve_units_dict["CAISO_reg_down"])

# define reserve service collection
reserveUp_services = collect(get_components(VariableReserve{ReserveUp}, sys))
reserveDown_services = collect(get_components(VariableReserve{ReserveDown}, sys))

# check to make sure reserve services were added
show_components(VariableReserve{ReserveUp}, sys)
show_components(VariableReserve{ReserveDown}, sys)
get_components(VariableReserve{ReserveUp}, sys) # retrieves an iterator of the reserve up services
get_components(VariableReserve{ReserveDown}, sys) # retrieves an iterator of the reserve down services

# spot check to make sure reserve memberships were assigned
# thermal standard
active_object = get_component(ThermalStandard, sys, "CAISO_CCGT1_PGE")
get_services(active_object)

# hydro
active_object = get_component(HydroDispatch, sys, "CAISO_Hydro_PGE")
get_services(active_object)

active_service = reserveUp_services[1]
for as_units in get_contributing_devices(sys, active_service)
    println(get_name(as_units))
end

active_service = reserveDown_services[1]
for as_units in get_contributing_devices(sys, active_service)
    println(get_name(as_units))
end

# now let's create PSY timeseries for regulation service
reg_ts_container_up, reg_ts_container_down, reg_reserve_df = create_reg_reserve_PSY_timeseries(demand_ts_df, gen_variability_ts_df, reserves_df, oprsv_zones, zone_dict, Renew_D_generators);

# Write regulation reserve components to CSV for visual inspection
CSV.write(joinpath(paths[:data_dir], "reg_reserve_components.csv"), reg_reserve_df);

#spot check the time series
active_ts = reg_ts_container_up["1998"]
active_ts.name
active_ts.data

active_ts = reg_ts_container_down["1998"]
active_ts.name
active_ts.data

# now that have our PSY timeseries, let's assign them to our reg service objects defined in the system
# ReserveUp
##########################
active_service = get_component(VariableReserve{ReserveUp}, sys, "CAISO_reg_up")
for (year, ts) in reg_ts_container_up
    # Only add up timeseries to up reserve service
    add_time_series!(sys, active_service, ts)
end

#check to make sure time series was added to requirement
show_time_series(active_service)
get_time_series_keys(active_service)
get_time_series_keys(active_service).ref
get_time_series_keys(active_service).size
get_time_series_array(SingleTimeSeries, active_service, "requirement_up_1998"; ignore_scaling_factors = true) # service was defined in natural units
get_time_series_array(SingleTimeSeries, active_service, "requirement_up_1998"; ignore_scaling_factors = false) # service was defined in natura units
get_requirement(active_service)
ts_active = get_time_series_array(SingleTimeSeries, active_service, "requirement_up_1999"; ignore_scaling_factors = true)
#set_requirement!(active_service, 5) # this jumps update the requirement attribute (has to be a constant)
get_requirement(active_service)

# ReserveDown
##########################
active_service = get_component(VariableReserve{ReserveDown}, sys, "CAISO_reg_down")
for (year, ts) in reg_ts_container_down
    # Only add up timeseries to up reserve service
    add_time_series!(sys, active_service, ts)
end

#check to make sure time series was added to requirement
show_time_series(active_service)
get_time_series_array(SingleTimeSeries, active_service, "requirement_down_1998"; ignore_scaling_factors = true) # service was defined in natural units
get_time_series_array(SingleTimeSeries, active_service, "requirement_down_1998"; ignore_scaling_factors = false) # service was defined in natura units


#= # remove time series from ReserveUp service objects   
for service in reserveUp_services
    # Get all time series keys for this service
    ts_keys = get_time_series_keys(service)
    # Loop through each time series key
    for ts_key in ts_keys
        ts_name = get_name(ts_key)
        # remove the time series
        remove_time_series!(sys, SingleTimeSeries, service, ts_name)
    end
end

# remove time series from ReserveDown service objects   
for service in reserveDown_services
    # Get all time series keys for this service
    ts_keys = get_time_series_keys(service)
    # Loop through each time series key
    for ts_key in ts_keys
        ts_name = get_name(ts_key)
        # remove the time series
        remove_time_series!(sys, SingleTimeSeries, service, ts_name)
    end
end =#    

##########################
# Define Fuel Time Series Info
##########################
# Process fuel price data
fuels_df = CSV.read(joinpath(paths[:data_dir], "system", "Fuels_data.csv"), DataFrame)
fuel_ts_df = process_fuel_data(fuels_df)

#= # Troubleshooting: write the fuels_df to a csv for visual inspection
CSV.write(joinpath(paths[:data_dir], "fuels_data.csv"), fuel_ts_df) =#

# Create fuel price timeseries by fuel type
fuel_price_ts_dict = create_fuel_price_PSY_timeseries(fuel_ts_df)

# spot check the fuel price timeseries
active_ts = fuel_price_ts_dict["CA_Natural_Gas"]
active_ts.name
active_ts.data

# now let's apply the fuel price forecast to each thermal generator based on resource type and region
for thermal_standard in ThermalStandard_generators
    
    # Get the fuel type for this generator
    fuel_type = string(get_fuel(thermal_standard))
    
    # For COAL and NATURAL_GAS, we need to consider regional differences
    if fuel_type in ["COAL", "NATURAL_GAS"]
        # Get the bus name and look up the region
        bus_name = get_name(get_bus(thermal_standard))
        region = bus_region_ba_df[bus_region_ba_df.bus .== bus_name, :].region[1]
        
        if isempty(region)
            @warn "No region found for bus $bus_name for generator $(get_name(thermal_standard))"
            continue
        end
        
        # Look up the corresponding fuel price time series based on region and fuel type
        matching_rows = fuel_mapping_df[
            (fuel_mapping_df.region .== region) .& 
            (fuel_mapping_df.sienna_fuel_type .== fuel_type), 
            :]
    else
        # For all other fuel types, we don't need to consider region
        matching_rows = fuel_mapping_df[
            fuel_mapping_df.sienna_fuel_type .== fuel_type, 
            :]
    end
    
    if nrow(matching_rows) == 0
        @warn "No fuel price mapping found for fuel type $fuel_type for generator $(get_name(thermal_standard))"
        continue
    end
    
    fuel_price_key = matching_rows[1, :fuel_price_fx]
    
    # Get the timeseries for this fuel type
    if haskey(fuel_price_ts_dict, fuel_price_key)
        ts = fuel_price_ts_dict[fuel_price_key]
        # Add the timeseries to the generator
        add_time_series!(sys, thermal_standard, ts)
    else
        @warn "No fuel price timeseries found for fuel price key $fuel_price_key (fuel type: $fuel_type) for generator $(get_name(thermal_standard))"
    end
end

# let's check our work
active_object = get_component(ThermalStandard, sys, "CAISO_CCGT1_PGE")
show_time_series(active_object) # should now see a fuel_price time series
get_name.(get_time_series_keys(active_object)) # retrieve the names of all defined time series 
get_time_series_array(SingleTimeSeries, active_object, "fuel_price"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_object, "fuel_price"; ignore_scaling_factors = false) # this should be the same as the first one (natural units)
active_object.operation_cost # we still see the default fuel cost (because haven't applied the fix yet)

# so now let's make the fuel price fx serve as the source of fuel costs for the PSI simulation(s)
update_TS_fuel_price!(sys, ThermalStandard_generators)

# you guessed it: let's check our work!
get_time_series(SingleTimeSeries, active_object, "fuel_price")
get_time_series_array(SingleTimeSeries, active_object, "fuel_price"; ignore_scaling_factors = true)
get_time_series_array(SingleTimeSeries, active_object, "fuel_price"; ignore_scaling_factors = false) # this should be equal to the previous cmd
# get_time_series_array(DeterministicSingleTimeSeries, active_object, "fuel_price")
active_object.operation_cost # should now see a pointer to the fuel_price timeseries as the defintion for ful_cost attribute

#= # remove all fuel_price_fx time series from system
for gen in get_components(x -> has_time_series(x), ThermalStandard, sys)
    remove_time_series!(sys, SingleTimeSeries, gen, "fuel_price")
end
 =#

##########################
# Define Outage Data (#TO-DO))
##########################
# SiennaPRASInterface will handle this

##########################
# Define PowerSimulations.jl (PSI) template and model 
##########################
# define run_type
run_type = "Monte_Carlo"

# define output path for simulation file
simulation_file_path = paths[:sienna_simulation_dir]

# Define the range of weather years based on run type
if run_type == "Deterministic"
    weather_years = [1998]  # Single year for deterministic
    sim_file_path = joinpath(paths[:sienna_simulation_dir], "deterministic")
elseif run_type == "Monte_Carlo" 
    weather_years = 2000:2001  # Range of years for Monte Carlo
    sim_file_path = joinpath(paths[:sienna_simulation_dir], "stochastic")
else
    @warn "Incorrect setting for run_type; $run_type is not a valid option"
    return
end

# Delete directory if it exists and create a fresh one
if ispath(sim_file_path)
    rm(sim_file_path, recursive=true, force=true)
end
mkpath(sim_file_path)

# Initialize dictionary to store simulation results
sim_results_dict = Dict{String, SimulationResults}()

# Create an empty model reference
template_uc = ProblemTemplate()

# Define non-weather related Device Models
###########################
# Storage
define_storage_model(template_uc)

# PHS
define_PHS_model(template_uc)

# Define branch model
###########################
define_branch_model(template_uc)

# Define network model
###########################
# CopperPlate
# define_CopperPlate_model(template_uc)

# AreaInterchange
define_AreaNetwork_model(template_uc)

# Loop through each weather year
for wy in weather_years

    # Remove all forecasts from the system
    remove_time_series!(sys, DeterministicSingleTimeSeries)
    
    # assign our generic "requirement" timeseries for our reserves 
    create_generic_requirement_reserveUp_timeseries(sys, wy);
    create_generic_requirement_reserveDown_timeseries(sys, wy);

    # assign our daily generic "hydro_budget" timeseries for our hydro units
    create_generic_daily_hydrobudget_timeseries(sys, wy, HydroDispatch_generators)

    # now let's create our timeseries forecasts
    transform_single_time_series!(sys, Hour(48), Hour(24))

#=     # Verify timeseries setup
    active_object = get_component(ThermalStandard, sys, "CAISO_CCGT1_PGE")
    ts_test = get_time_series(DeterministicSingleTimeSeries, active_object, "fuel_price")
    horizon = get_horizon(ts_test)
    interval = get_interval(ts_test) =#

    # Define weather-dependent Device Models
    define_load_model(template_uc, wy)
    define_thermal_model(template_uc, wy)
    define_hydro_model(template_uc, wy)
    define_renewable_dispatch_model(template_uc, wy)
    define_renewable_non_dispatch_model(template_uc, wy)

    # define the service models
    ###########################
    #define_RegUp_service_model(template_uc) # remember: we already updated our timeseries for the active WY
    #define_RegDown_service_model(template_uc)

    # Assign simulation name
    decision_name = "$(lowercase(run_type))_$wy"

    # Initialize decision model
    UC_decision = DecisionModel(
        template_uc,
        sys;
        name = decision_name,
        optimizer = optimizer_with_attributes(Gurobi.Optimizer, "MIPGap" => 1e-2),
        system_to_file = false,
        initialize_model = true,
        optimizer_solve_log_print = true,
        direct_mode_optimizer = true,
        rebuild_model = false,
        store_variable_names = true,
        calculate_conflict = true,
        export_optimization_model = false,
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
        name = "sim_$wy",
        steps = 3,  # Steps in your simulation
        models = sim_model,
        sequence = sim_sequence,
        simulation_folder = sim_file_path,
    )

    # simulation_folder = mktempdir(paths[:sienna_simulation_dir], cleanup = true),

    # Build and execute simulation
    build!(sim; console_level = Logging.Info)
    execute!(sim, enable_progress_bar = true)

    # Store simulation results in dictionary
    sim_results = SimulationResults(sim)
    sim_results_dict[decision_name] = sim_results

    # Print completion message
    println("Simulation completed for: $decision_name")
end

###########################
# Export the Results
###########################
# Separate loop for processing and exporting results
for (decision_name, sim_results) in sim_results_dict
    println("\nProcessing results for: $decision_name")

    #retrieve the weather year from the decision name
    wy = parse(Int, decision_name[end-3:end])

    #define the results file path based on the run type
    if run_type == "Deterministic"
        # define output path for results file
        results_file_path = joinpath(paths[:PSI_results_dir], "results_$wy")
    elseif run_type == "Monte_Carlo" 
        results_file_path = joinpath(paths[:PRAS_results_dir], "results_$wy")
    else
        @warn "Incorrect setting for run_type; $run_type is not a valid option"
    end
    
    # Check if the results directory exists, if not create it
    if !ispath(results_file_path)
        mkpath(results_file_path)
    end
    
    # Get results for this weather year
    results = get_decision_problem_results(sim_results, decision_name)
    
    # Input TimeSeries Parameters
    load_parameter = read_realized_parameter(results, "ActivePowerTimeSeriesParameter__PowerLoad")
    thermal_parameter = read_realized_parameter(results, "ActivePowerTimeSeriesParameter__ThermalStandard")
    renewDispatch_parameter = read_realized_parameter(results, "ActivePowerTimeSeriesParameter__RenewableDispatch")
    renewNonDispatch_parameter = read_realized_parameter(results, "ActivePowerTimeSeriesParameter__RenewableNonDispatch")
    hydro_dispatch_parameter = read_realized_parameter(results, "ActivePowerTimeSeriesParameter__HydroDispatch")
    hydro_budget_parameter = read_realized_parameter(results, "EnergyBudgetTimeSeriesParameter__HydroDispatch")
    PHS_output_parameter = read_realized_parameter(results, "OutflowTimeSeriesParameter__HydroPumpedStorage")
    PHS_input_parameter = read_realized_parameter(results, "InflowTimeSeriesParameter__HydroPumpedStorage")

    # Realized Variables
    thermal_active_power = read_realized_variable(results, "ActivePowerVariable__ThermalStandard")
    renewDispatch_active_power = read_realized_variable(results, "ActivePowerVariable__RenewableDispatch")
    battery_charge = read_realized_variable(results, "ActivePowerInVariable__EnergyReservoirStorage")
    battery_discharge = read_realized_variable(results, "ActivePowerOutVariable__EnergyReservoirStorage")
    battery_energy = read_realized_variable(results, "EnergyVariable__EnergyReservoirStorage")
    hydro_active_power = read_realized_variable(results, "ActivePowerVariable__HydroDispatch")
    PHS_charge = read_realized_variable(results, "ActivePowerInVariable__HydroPumpedStorage")
    PHS_discharge = read_realized_variable(results, "ActivePowerOutVariable__HydroPumpedStorage")
    PHS_spillage = read_realized_variable(results, "WaterSpillageVariable__HydroPumpedStorage")
    PHS_UpperReservoir = read_realized_variable(results, "HydroEnergyVariableUp__HydroPumpedStorage")
    PHS_LowerReservoir = read_realized_variable(results, "HydroEnergyVariableDown__HydroPumpedStorage")
    PHS_reservation = read_realized_variable(results, "ReservationVariable__HydroPumpedStorage")

    # Combine results
    gen_power = hcat(thermal_active_power, select(renewDispatch_active_power, Not(1)), select(hydro_active_power, Not(1)))
    storage_discharge_power = hcat(PHS_discharge, select(battery_discharge, Not(1)))
    storage_charge_power = hcat(PHS_charge, select(battery_charge, Not(1)))
    AreaInterchange_flow = read_realized_variable(results, "FlowActivePowerVariable__AreaInterchange")
    power_balance = read_realized_expression(results, "ActivePowerBalance__Area")

    # Production Costs
    pc_thermal = read_realized_expression(results, "ProductionCostExpression__ThermalStandard")
    pc_renewable = read_realized_expression(results, "ProductionCostExpression__RenewableDispatch")
    pc_hydro = read_realized_expression(results, "ProductionCostExpression__HydroDispatch")
    pc_PHS = read_realized_expression(results, "ProductionCostExpression__HydroPumpedStorage")
    pc_all = hcat(pc_thermal, select(pc_renewable, Not(1)), select(pc_hydro, Not(1)), select(pc_PHS, Not(1)))
    fuel_consumption_thermal = read_realized_expression(results, "FuelConsumptionExpression__ThermalStandard")

    # Auxiliary variables
    Energy_PHS = read_realized_variable(results, "HydroEnergyOutput__HydroPumpedStorage")
    Energy_Hydro = read_realized_variable(results, "HydroEnergyOutput__HydroDispatch")
    Energy_Battery = read_realized_variable(results, "StorageEnergyOutput__EnergyReservoirStorage")

    println("Exporting results to CSV files...")

    # Export results to CSV files
    CSV.write(joinpath(results_file_path, "power_load_parameters.csv"), load_parameter)
    CSV.write(joinpath(results_file_path, "thermal_parameters.csv"), thermal_parameter)
    CSV.write(joinpath(results_file_path, "FTM_renewable_parameters.csv"), renewDispatch_parameter)
    CSV.write(joinpath(results_file_path, "BTM_renewable_parameters.csv"), renewNonDispatch_parameter)
    CSV.write(joinpath(results_file_path, "hydro_parameter.csv"), hydro_dispatch_parameter)
    CSV.write(joinpath(results_file_path, "hydro_budget_parameter.csv"), hydro_budget_parameter)
    CSV.write(joinpath(results_file_path, "PHS_output_parameter.csv"), PHS_output_parameter)
    CSV.write(joinpath(results_file_path, "PHS_input_parameter.csv"), PHS_input_parameter)

    CSV.write(joinpath(results_file_path, "FTM_generator_power.csv"), gen_power)
    CSV.write(joinpath(results_file_path, "storage_charge.csv"), storage_charge_power)
    CSV.write(joinpath(results_file_path, "storage_discharge.csv"), storage_discharge_power)
    CSV.write(joinpath(results_file_path, "AreaInterchange_flow.csv"), AreaInterchange_flow)
    CSV.write(joinpath(results_file_path, "power_balance.csv"), power_balance)
    CSV.write(joinpath(results_file_path, "production_costs.csv"), pc_all)
    CSV.write(joinpath(results_file_path, "fuel_consumption_thermal.csv"), fuel_consumption_thermal)

    # PHS Specific
    CSV.write(joinpath(results_file_path, "PHS_UpperReservoir.csv"), PHS_UpperReservoir)
    CSV.write(joinpath(results_file_path, "PHS_LowerReservoir.csv"), PHS_LowerReservoir)
    CSV.write(joinpath(results_file_path, "PHS_spillage.csv"), PHS_spillage)
    CSV.write(joinpath(results_file_path, "PHS_reservation.csv"), PHS_reservation)
    CSV.write(joinpath(results_file_path, "Energy_PHS.csv"), Energy_PHS)

    # Hydro Specific
    CSV.write(joinpath(results_file_path, "Energy_Hydro.csv"), Energy_Hydro)

    # Battery Storage Specific
    CSV.write(joinpath(results_file_path, "Energy_Battery.csv"), Energy_Battery)

end

end # module
