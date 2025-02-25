module GenX_to_Sienna

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

# parameters

data_directory = "/scratch/gpfs/gm1710/SCP_GenX_Case_FullWY" 
output_directory_base = "/scratch/gpfs/gm1710/SCP_Sienna_Outputs/results_"
results_folder_name = "results_resolverepperiods"

global start_time::Int
start_time=1
for i in 0:22
    global start_time
    year=1998+i
    output_directory = output_directory_base*string(year)
    mkdir(output_directory)
    if (mod(year,4)==0) && (year!=2020)
        end_time = start_time+8784-1
    else
        end_time = start_time+8760-1
    end
    num_hours = end_time-start_time+1
    timestamps = range(DateTime(string(year)*"-01-01T00:00:00"), step = Dates.Hour(1), length = num_hours)
    horizon=num_hours

    # read in relevant CSV files

    if isfile(joinpath(data_directory, "results", "network_expansion.csv"))
        network_expansion=true
    else
        network_expansion=false
    end

    if isfile(joinpath(data_directory, "resources", "Must_run.csv"))
        must_run=true
    else
        must_run=false
    end

    storage_df = CSV.read(joinpath(data_directory, "resources", "Storage.csv"), DataFrame) ;
    thermal_df = CSV.read(joinpath(data_directory, "resources", "Thermal.csv"), DataFrame) ;
    vre_df = CSV.read(joinpath(data_directory, "resources", "Vre.csv"), DataFrame) ;
    hydro_df = CSV.read(joinpath(data_directory, "resources", "Hydro.csv"), DataFrame) ;
    capacity_df = CSV.read(joinpath(data_directory, results_folder_name, "capacity.csv"), DataFrame);
    demand_df = CSV.read(joinpath(data_directory, "system", "Demand_data.csv"), DataFrame) ;
    network_df = CSV.read(joinpath(data_directory, "system", "Network.csv"), DataFrame) ;
    reserves_df = CSV.read(joinpath(data_directory, "system", "Operational_reserves.csv"), DataFrame) ;
    if network_expansion
        network_expansion_df = CSV.read(joinpath(data_directory, "results", "network_expansion.csv"), DataFrame) ;
    end
    if must_run
        must_run_df = CSV.read(joinpath(data_directory, "resources", "Must_run.csv"), DataFrame) ;
    end
    fuels_df = CSV.read(joinpath(data_directory, "system", "Fuels_data.csv"), DataFrame) ;
    gen_variability_df = CSV.read(joinpath(data_directory, "system", "Generators_variability.csv"), DataFrame) ;
    energy_budget_df = CSV.read(joinpath(data_directory, "system", "Hourly_energy_budget.csv"), DataFrame) ;
    mt_df = CSV.read(joinpath(data_directory, "MoverTypesMapping.csv"), DataFrame) ;
    fm_df = CSV.read(joinpath(data_directory, "FuelMapping.csv"), DataFrame) ;
    sm_df = CSV.read(joinpath(data_directory, "StorageMapping.csv"), DataFrame) ;
    rm_df = CSV.read(joinpath(data_directory, "RenewableMapping.csv"), DataFrame) ;

    # create 3 dictionaries, one for fuel types, one for mover types, one for storage types
    mover_dict = Dict((row.Key) => row.Value for row in eachrow(mt_df)) ;
    fuel_dict = Dict((row.Key) => row.Value for row in eachrow(fm_df)) ;
    storage_dict = Dict((row.Key) => row.Value for row in eachrow(sm_df)) ;
    renewable_dict = Dict((row.Key) => row.Value for row in eachrow(rm_df)) ;

    #fuel costs:
    column_names = names(fuels_df)
    columns_to_read = column_names[2:end]
    fuel_prices = Dict{String, Union{Float64, Missing}}()
    for col in columns_to_read
        average = mean(fuels_df[2:end, Symbol(col)])  # Ignore the first row
        fuel_prices[String(col)] = average
    end


    # adding nodes
    nodes=[]
    for i in 1:count(!ismissing, network_df[:, "Network_zones"])
        bus=ACBus(;
            number=i,
            name= network_df[i, :1],
            bustype= i == 1 ? "REF" : "PQ",
            angle=0,
            magnitude=0,
            voltage_limits=nothing,
            base_voltage=230
        )
        push!(nodes,bus)
    end
    nodes=Vector{ACBus}(nodes)

    # make system
    sys_base_power=100.0
    sys=System(sys_base_power,nodes)
    set_units_base_system!(sys, "NATURAL_UNITS") ;

    areas=[]
    for i in 1:count(!ismissing, network_df[:, "Network_zones"])
        area_factor=Area(;name=network_df[i, :1])
        push!(areas,area_factor)
        set_area!(nodes[i], areas[i])
    end
    add_components!(sys,areas)

    # create line vector
    lines=[]
    for i in 1:count(!ismissing, network_df[:, "Network_Lines"])
        existing_cap = network_df[i, :Line_Max_Flow_MW]
        if network_expansion
            new_cap = network_expansion_df[network_expansion_df.Line .== network_df[i, :Network_Lines], "New_Trans_Capacity"][1]
        else
            new_cap=0.0
        end
        line = Line(;
            name = string(network_df[i, :Network_Lines]),
            available = true,
            active_power_flow = 0.0,
            reactive_power_flow = 0.0,
            arc=Arc(from = nodes[network_df[i, :Start_Zone]], to = nodes[network_df[i, :End_Zone]]),
            r=0.0,
            x=0.0,
            b= (from = 0.0, to = 0.0),
            rating = (existing_cap+new_cap)/sys_base_power,
            angle_limits = (0.0, 0.0),
        )
        push!(lines,line)  
    end
    add_components!(sys,lines) # add lines vector to system

    ais=[]
    for i in 1:count(!ismissing, network_df[:,"Network_Lines"])
        ai=AreaInterchange(;
            name = string(network_df[i, :Network_Lines]),
            available = true,
            active_power_flow = 0.0,
            flow_limits=(from_to=get_rating(lines[i])*network_df[i, :Profile_Forward]/sys_base_power,to_from=get_rating(lines[i])*network_df[i, :Profile_Reverse]/sys_base_power),
            from_area=areas[network_df[i, :Start_Zone]],
            to_area=areas[network_df[i, :End_Zone]],
        )
        push!(ais,ai)
    end
    add_components!(sys,ais)

    # initiate list of resources contributing to frequency regulation and operating reserves
    vector_of_regulation_resources = []
    vector_of_reserves_resources = []

    # thermal generators
    thermal_gen=[]
    for i in 1:count(!ismissing, thermal_df[:, "Resource"])
        resource_name = thermal_df[i, :Resource]
        capacity_mw = capacity_df[capacity_df.Resource .== resource_name, :][1, :EndCap]
        if "Heat_Rate_MMBTU_per_MWh" in names(thermal_df)
            heat_rate_curve = LinearCurve(thermal_df[i, :Heat_Rate_MMBTU_per_MWh],0.0)
        else
            heat_rate_curve = LinearCurve(thermal_df[i, :PWFU_Heat_Rate_MMBTU_per_MWh_1],thermal_df[i, :PWFU_Fuel_Usage_Zero_Load_MMBTU_per_h])
        end
        if "PWFU_Heat_Rate_MMBTU_per_MWh_2" in names(thermal_df)
            throw(ErrorException("More than one PWFU segment not supported."))
        end
        thermal = ThermalStandard(;
            name = resource_name,
            available = true,
            status = true,
            bus = nodes[thermal_df[i, :Zone]],
            active_power = 0,
            reactive_power = 0,
            rating = 1,
            active_power_limits = (min = thermal_df[i, :Min_Power], max = 1),
            reactive_power_limits = nothing,
            ramp_limits = (up=thermal_df[i, :Ramp_Up_Percentage]/60, down = thermal_df[i, :Ramp_Dn_Percentage]/60),
            operation_cost = ThermalGenerationCost(
                variable = FuelCurve(
                    value_curve = heat_rate_curve,
                    fuel_cost = fuel_prices[string(thermal_df[i, :Fuel])]
                    ),
                fixed = 0,
                start_up = thermal_df[i, :Start_Cost_per_MW]*capacity_mw,
                shut_down = 0.0,
                ),
            base_power = capacity_mw,
            time_limits = (up = thermal_df[i, :Up_Time], down = thermal_df[i, :Down_Time]),
            must_run = false,
            prime_mover_type = getproperty(PrimeMovers,Symbol(mover_dict[thermal_df[i, :Resource]])),
            fuel = fuel_dict[string(thermal_df[i, :Fuel])],
        )
        push!(thermal_gen, thermal)
        if thermal_df[i,:Reg_Max]>0
            push!(vector_of_regulation_resources, thermal)
        end
        if thermal_df[i,:Rsv_Max]>0
            push!(vector_of_reserves_resources, thermal)
        end
    end
    add_components!(sys,thermal_gen) 

    # renewable generators
    for i in 1:count(!ismissing, vre_df[:, "Resource"])
        resource_name = vre_df[i, :Resource]
        renewable = RenewableDispatch(;
            name = resource_name,
            available = true,
            bus = nodes[vre_df[i, :Zone]],
            active_power = 0.0,
            reactive_power = 0.0,
            rating = 1,
            prime_mover_type = getproperty(PrimeMovers,Symbol(mover_dict[resource_name])),
            reactive_power_limits = (min = 0.0, max = 0.0),
            power_factor = 1.0,
            operation_cost = RenewableGenerationCost(CostCurve(LinearCurve(vre_df[i, :Var_OM_Cost_per_MWh]))),
            base_power = capacity_df[capacity_df.Resource .== resource_name, :][1, :EndCap]
            )
        add_component!(sys,renewable)
        if vre_df[i,:Reg_Max]>0
            push!(vector_of_regulation_resources, renewable)
        end
        if vre_df[i,:Rsv_Max]>0
            push!(vector_of_reserves_resources, renewable)
        end
        ren_values = gen_variability_df[start_time:end_time, Symbol(resource_name)]
        ren_timearray = TimeArray(timestamps, ren_values);
        ren_time_series = SingleTimeSeries(
            name = "max_active_power",
            data = ren_timearray;
            scaling_factor_multiplier = get_max_active_power
        );
        add_time_series!(sys, renewable, ren_time_series);
    end

    # must run generators (assume only Customer PV)
    for i in 1:count(!ismissing, must_run_df[:, "Resource"])
        resource_name = must_run_df[i, :Resource]
        if resource_name!="Customer_PV"
            throw(ErrorException("Only Customer PV is supported as a must-run generator."))
        end
        renewable = RenewableNonDispatch(;
            name = resource_name,
            available = true,
            bus = nodes[must_run_df[i, :Zone]],
            active_power = 0.0,
            reactive_power = 0.0,
            rating = 1.0,
            prime_mover_type = getproperty(PrimeMovers,Symbol("PVe")),
            power_factor = 1.0,
            base_power = capacity_df[capacity_df.Resource .== resource_name, :][1, :EndCap]
            )
        add_component!(sys,renewable)
        ren_values = gen_variability_df[start_time:end_time, Symbol(resource_name)]
        ren_timearray = TimeArray(timestamps, ren_values);
        ren_time_series = SingleTimeSeries(
            name = "max_active_power",
            data = ren_timearray;
            scaling_factor_multiplier = get_max_active_power
        );
        add_time_series!(sys, renewable, ren_time_series);
    end

    # hydro
    for i in 1:count(!ismissing, hydro_df[:, "Resource"])
        resource_name = hydro_df[i, :Resource]
        hydro = ThermalStandard(;
            name = resource_name,
            available = true,
            status=true,
            bus = nodes[hydro_df[i, :Zone]],
            active_power = 0.0,
            reactive_power = 0.0,
            rating = 1.0,
            prime_mover_type = getproperty(PrimeMovers,Symbol("CC")),
            active_power_limits = (min = hydro_df[i, :Min_Power], max = 1),
            reactive_power_limits = (min = 0.0, max = 0.0),
            ramp_limits = (up=hydro_df[i, :Ramp_Up_Percentage]/60, down = hydro_df[i, :Ramp_Dn_Percentage]/60),
            time_limits = nothing,
            base_power = capacity_df[capacity_df.Resource .== resource_name, :][1, :EndCap],
            operation_cost = ThermalGenerationCost(
            variable=CostCurve(LinearCurve(hydro_df[i, :Var_OM_Cost_per_MWh])),
            fixed=0.0,
            start_up=0.0,
            shut_down=0.0),
        )
        add_component!(sys,hydro)
        if hydro_df[i,:Reg_Max]>0
            push!(vector_of_regulation_resources, hydro)
        end
        if hydro_df[i,:Rsv_Max]>0
            push!(vector_of_reserves_resources, hydro)
        end
        ren_values = gen_variability_df[start_time:end_time, Symbol(resource_name)]
        energy_budget = energy_budget_df[start_time:end_time, Symbol(resource_name)]
        cf_to_assign = min(ren_values,energy_budget)
        ren_timearray = TimeArray(timestamps, cf_to_assign);
        ren_time_series = SingleTimeSeries(
            name = "max_active_power",
            data = ren_timearray;
            scaling_factor_multiplier = get_max_active_power
        );
        add_time_series!(sys, hydro, ren_time_series);
    end

    # storage
    for i in 1:count(!ismissing, storage_df[:, "Resource"])
        resource_name = storage_df[i, :Resource]
        storage_device=EnergyReservoirStorage(;
            name = resource_name,
            available = true,
            bus = nodes[storage_df[i, :Zone]],
            prime_mover_type = getproperty(PrimeMovers,Symbol(mover_dict[resource_name])),
            storage_technology_type = getproperty(StorageTech,Symbol(storage_dict[resource_name])),
            storage_capacity = capacity_df[capacity_df.Resource .== resource_name, :][1, :EndEnergyCap],
            storage_level_limits = (min = 0.0, max = 1.0),
            initial_storage_capacity_level = 0.0,
            rating = 1.0,
            active_power = 0.0,
            input_active_power_limits = (min = 0.0, max = 1.0),
            output_active_power_limits = (min = 0.0, max = 1.0),
            efficiency = (in = storage_df[i, :Eff_Up], out = storage_df[i, :Eff_Down]),
            reactive_power = 0,
            reactive_power_limits = (min = 0.0, max = 0.0),
            base_power = capacity_df[capacity_df.Resource .== resource_name, :][1, :EndCap],
            operation_cost = StorageCost(charge_variable_cost = CostCurve(LinearCurve(storage_df[i, :Var_OM_Cost_per_MWh_In])), 
                discharge_variable_cost = CostCurve(LinearCurve(storage_df[i, :Var_OM_Cost_per_MWh])), 
                fixed = 0.0, start_up = 0.0, shut_down = 0.0, energy_shortage_cost = 0.0, 
                energy_surplus_cost = 0.0),
            )
        add_component!(sys, storage_device)
        if storage_df[i,:Reg_Max]>0
            push!(vector_of_regulation_resources, storage_device)
        end
        if storage_df[i,:Rsv_Max]>0
            push!(vector_of_reserves_resources, storage_device)
        end
    end

    # add NSE as a thermal generator for each zone

    VoLL = demand_df[1,"Voll"]
    for i in 1:count(!ismissing, network_df[:, "Network_zones"])
        nse = ThermalStandard(;
            name = "VoLL_"*string(i),
            available = true,
            status = true,
            bus = nodes[i],
            active_power = 0,
            reactive_power = 0,
            rating = 1,
            active_power_limits = (min = 0, max = 1),
            reactive_power_limits = nothing,
            ramp_limits = (up=1000, down = 1000),
            operation_cost=ThermalGenerationCost(
                variable=CostCurve(LinearCurve(VoLL)),
                fixed=0.0,
                start_up=0.0,
                shut_down=0.0
            ),
            base_power = 1e6,
            must_run = false,
            prime_mover_type = getproperty(PrimeMovers,Symbol("CC")),
        )
        add_component!(sys, nse)
    end

    # adding load to each zone:
    for i in 1:count(!ismissing, network_df[:, "Network_zones"])
        zone_name = "Load_z"*string(i)
        power_load=PowerLoad(;
            name=zone_name,
            available=true,
            bus=nodes[i],
            active_power=1e9,
            reactive_power=1e9,
            base_power=1,
            max_active_power=1,
            max_reactive_power=0
            )
        add_component!(sys, power_load)
        data = TimeArray(timestamps, Float64.(demand_df[start_time:end_time,Symbol("Demand_MW_z" * "" * string(i))]))
        time_series = SingleTimeSeries(name="max_active_power", data=data,scaling_factor_multiplier=get_max_active_power)
        add_time_series!(sys,get_component(PowerLoad, sys, zone_name),time_series)
    end

    # add frequency regulation reserves to relevant generators

    regulation = VariableReserve{ReserveUp}(
        name="regulation",
        available=true,
        time_frame=0.0,
        requirement=1/sys_base_power # 1 MW, will later be scaled
        )
    add_service!(sys,regulation,vector_of_regulation_resources)

    # add reserve requirement based on demand
    reg_requirement=reserves_df[1,:Reg_Req_Percent_Demand]
    reg_zone=reserves_df[1,:OpRsv_Zone]
    demand_timeseries = Float64.(demand_df[start_time:end_time,Symbol("Demand_MW_z" * "" * string(reg_zone))])
    reserve_requirement = demand_timeseries * reg_requirement
    reserve_timearray = TimeArray(timestamps, reserve_requirement);
    reserve_timeseries = SingleTimeSeries(;
            name = "requirement",
            data = reserve_timearray,
            scaling_factor_multiplier = get_requirement,
        );
    add_time_series!(sys, regulation, reserve_timeseries);

    transform_single_time_series!(sys,Dates.Hour(horizon),Dates.Hour(horizon)) 

    template_uc=ProblemTemplate()
    set_device_model!(template_uc,ThermalStandard,ThermalDispatchNoMin)
    set_device_model!(template_uc,PowerLoad,StaticPowerLoad)
    set_device_model!(template_uc, RenewableDispatch, RenewableFullDispatch)
    set_device_model!(template_uc, RenewableNonDispatch, FixedOutput)
    set_device_model!(template_uc, Line, StaticBranch)
    set_device_model!(template_uc,DeviceModel(
            EnergyReservoirStorage,
            StorageDispatchWithReserves;
            attributes=Dict("reservation" => true,
                "cycling_limits" => false,
                "energy_target" => false,
                "complete_coverage" => false,
                "regularization" => true,
            ),
            use_slacks=false,
        )
    )
    set_network_model!(template_uc, NetworkModel(AreaBalancePowerModel))
    set_device_model!(template_uc, AreaInterchange, StaticBranch)
    set_service_model!(template_uc, VariableReserve{ReserveUp}, RangeReserve)
    problem=DecisionModel(template_uc,sys;optimizer=optimizer_with_attributes(Gurobi.Optimizer),horizon=Dates.Hour(horizon))

    println("building model")
    build!(problem,output_dir=mktempdir())
    println("solving model")
    solve!(problem)
    println("getting outputs and writing")
    res=OptimizationProblemResults(problem) ;
    gen = get_generation_data(res) ;

    # get model data
    charge_energyreservoirstorage=gen.data[:ActivePowerInVariable__EnergyReservoirStorage]
    discharge_energyreservoirstorage=gen.data[:ActivePowerOutVariable__EnergyReservoirStorage]
    activepower_thermalstandard=gen.data[:ActivePowerVariable__ThermalStandard]
    df_ren_dis=gen.data[:ActivePowerVariable__RenewableDispatch]
    demand=demand_df[start_time:end_time,(end-(length(areas)-1)):end]
    curtailment=gen.data[:ActivePowerVariable__RenewableDispatch__Curtailment]
    flows=read_variable(res,"FlowActivePowerVariable__AreaInterchange")

    CSV.write(joinpath(output_directory,"charge_energyreservoirstorage.csv"),charge_energyreservoirstorage)
    CSV.write(joinpath(output_directory,"discharge_energyreservoirstorage.csv"),discharge_energyreservoirstorage)
    CSV.write(joinpath(output_directory,"activepower_thermalstandard.csv"),activepower_thermalstandard)
    CSV.write(joinpath(output_directory,"df_ren_dis.csv"),df_ren_dis)
    CSV.write(joinpath(output_directory,"demand.csv"),demand)
    CSV.write(joinpath(output_directory,"curtailment.csv"),curtailment)
    CSV.write(joinpath(output_directory,"flows.csv"),flows)

    start_time=end_time+1

end

end # module GenX_to_Sienna
