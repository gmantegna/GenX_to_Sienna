# GenX to Sienna Conversion Tool

This project converts GenX optimization model outputs into the Sienna ecosystem (PowerSystems.jl and PowerSimulations.jl) for advanced power system analysis and simulation.

## Overview

The **GenX to Sienna** conversion tool takes GenX capacity expansion planning results and transforms them into a comprehensive PowerSystems.jl system model, then executes production cost modeling simulations using PowerSimulations.jl. This enables detailed operational analysis of power system portfolios optimized through GenX.

## Project Structure

```
src/
├── GenX_to_Sienna.jl           # Main orchestration file
├── helpers.jl                  # Utility functions and data processing
├── system_builder.jl           # PowerSystems component creation
├── time_series_builder.jl      # Time series data management
└── simulation_builder.jl       # PowerSimulations model definitions
```

## Workflow Overview

### 1. Main Orchestration (`GenX_to_Sienna.jl`)

The main file coordinates the entire conversion process through the following stages:

#### **System Initialization**
- Sets up directory paths and logging using `helpers.jl`
- Creates a base PowerSystems.jl `System` object with 1.0 MVA base power
- Configures natural units (vs. per-unit) for easier data interpretation

#### **Network Topology Construction**
- **Buses**: Creates `ACBus` objects for each GenX zone using `system_builder.jl`
- **Areas**: Defines `Area` objects and assigns buses to areas
- **Lines**: Builds `Line` objects from GenX network topology data
- **Area Interchanges**: Creates `AreaInterchange` objects for inter-area power flows
- **Transmission Interfaces**: Defines `TransmissionInterface` objects for simultaneous flow constraints

#### **Load Definition**
- Creates `PowerLoad` objects for each area using demand data
- Processes demand time series through `helpers.jl` data processing functions

#### **Generation Portfolio**
- **Thermal Generators**: Creates `ThermalStandard` objects with operational parameters
- **Renewable Dispatch**: Creates `RenewableDispatch` objects for utility-scale VRE
- **Renewable Non-Dispatch**: Creates `RenewableNonDispatch` objects for behind-the-meter resources
- **Hydro**: Creates `HydroDispatch` objects with energy budget constraints
- **Battery Storage**: Creates `EnergyReservoirStorage` objects
- **Pumped Hydro Storage**: Creates `HydroPumpedStorage` objects

#### **Time Series Integration**
- Attaches operational time series to all components using `time_series_builder.jl`
- Handles multiple weather years for stochastic analysis
- Creates forecast windows for rolling horizon simulations

#### **Ancillary Services**
- Defines regulation reserve services (`VariableReserve`) for up and down reserves
- Assigns eligible resources to reserve markets

#### **Simulation Execution**
- Builds PowerSimulations.jl templates using `simulation_builder.jl`
- Executes simulations across multiple weather years
- Supports both deterministic and Monte Carlo analysis modes

#### **Results Export**
- Processes and exports simulation results to CSV files
- Organizes outputs by weather year and analysis type

### 2. Supporting Modules

#### **`helpers.jl` - Utility Functions**
- **Path Management**: `initialize_paths_and_inputs()` sets up directory structure
- **Data Processing**: Functions to process GenX CSV outputs into Julia DataFrames
  - `process_demand_data()` - Handles load time series
  - `process_fuel_data()` - Processes fuel price data
  - `process_generator_variability_data()` - Handles renewable profiles
  - `process_hydro_budget_data()` - Processes hydro energy budgets
- **Output Functions**: Creates parameter summary CSV files for validation

#### **`system_builder.jl` - PowerSystems Component Creation**
- **Network Components**:
  - `create_buses()` - Creates ACBus objects from GenX zones
  - `create_areas()` - Creates Area objects and bus assignments
  - `create_lines()` - Creates Line objects with thermal ratings
  - `create_area_interchanges()` - Creates AreaInterchange objects
  - `create_transmission_interfaces()` - Creates TransmissionInterface objects

- **Generation Components**:
  - `create_power_loads()` - Creates PowerLoad objects
  - `create_ThermalStandard_objects()` - Creates thermal generators with heat rates and operational costs
  - `create_VRE_objects()` - Creates renewable dispatch generators
  - `create_btm_objects()` - Creates behind-the-meter renewable resources
  - `create_Hydro_objects()` - Creates hydro generators with energy budgets
  - `create_storage_objects()` - Creates battery storage systems
  - `create_PHS_objects()` - Creates pumped hydro storage systems

- **Ancillary Services**:
  - `create_CAISO_reg_reserve_services()` - Creates reserve service objects
  - `create_CAISO_reg_reserve_units()` - Assigns eligible resources to reserves

#### **`time_series_builder.jl` - Time Series Management**
- **Load Time Series**: `create_demand_PSY_timeseries()` - Creates demand forecasts
- **Generation Time Series**:
  - `create_Renew_D_PSY_timeseries()` - Renewable dispatch availability profiles
  - `create_ThermalStandard_PSY_timeseries()` - Thermal availability profiles
  - `create_Hydro_PSY_timeseries()` - Hydro availability profiles
  - `create_Hydro_Budget_PSY_timeseries()` - Hydro energy budget time series
- **Reserve Time Series**: `create_reg_reserve_PSY_timeseries()` - Reserve requirement forecasts
- **Fuel Price Time Series**: `create_fuel_price_PSY_timeseries()` - Variable fuel costs
- **Forecast Management**: Functions to transform single time series into rolling forecasts

#### **`simulation_builder.jl` - PowerSimulations Model Definitions**
- **Device Models**:
  - `define_thermal_model()` - Unit commitment formulation for thermal generators
  - `define_renewable_dispatch_model()` - Economic dispatch for renewables
  - `define_hydro_model()` - Hydro dispatch with energy budgets
  - `define_storage_model()` - Storage optimization with cycling constraints
  - `define_PHS_model()` - Pumped hydro storage dispatch
- **Network Models**:
  - `define_AreaNetwork_model()` - Zonal network representation
  - `define_branch_model()` - Transmission line flow constraints
- **Service Models**: Reserve requirement and allocation models

## Key Features

### Multi-Year Analysis
- Supports both deterministic (single weather year) and stochastic (multiple weather years) analysis
- Handles leap year adjustments in time series data
- Maintains temporal consistency across different data sources

### Comprehensive Modeling
- **Full Generation Portfolio**: Thermal, renewable, hydro, and storage resources
- **Transmission Constraints**: Line limits, area interchange limits, and interface constraints
- **Ancillary Services**: Regulation reserves with dynamic requirements
- **Operational Details**: Unit commitment, ramping constraints, minimum generation levels

### Data Integration
- Seamless conversion from GenX CSV outputs to PowerSystems.jl components
- Automatic fuel price time series integration
- Comprehensive parameter validation and export capabilities

## Usage

1. **Setup**: Ensure GenX output files are in the expected directory structure
2. **Configuration**: Modify run type (Deterministic vs Monte Carlo) and weather year ranges
3. **Execution**: Run `GenX_to_Sienna.jl` to perform the full conversion and simulation
4. **Results**: Access processed results in organized CSV output files

## Dependencies

- **PowerSystems.jl**: Core power system data modeling
- **PowerSimulations.jl**: Production cost modeling and optimization
- **PowerAnalytics.jl**: Results processing and analysis
- **HydroPowerSimulations.jl**: Advanced hydro modeling capabilities
- **StorageSystemsSimulations.jl**: Enhanced storage modeling
- **Gurobi.jl**: Commercial optimization solver interface

## Output Structure

Results are organized by analysis type and weather year:
```
Sienna_Outputs/
├── output_results/
│   ├── PowerSimulations/    # Deterministic results
│   └── PRAS/               # Monte Carlo results
└── simulation_files/       # Simulation metadata and logs
```

This tool bridges the gap between long-term capacity planning (GenX) and detailed operational analysis (Sienna), enabling comprehensive power system studies that span from investment decisions to real-time operations. 