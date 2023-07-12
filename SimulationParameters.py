import CoolProp.CoolProp as CP
import pandas as pd
from PPIs import PPIs

class SimulationParameters:
    def __init__(self):
        # depth [m], positive
        self.depth = 6000

        # mass flowrate [kg/s], positive
        self.m_dot = None

        # inside well radius [m] ~18 3/8"
        self.well_radius = 0.5087

        # side stream inside radius [m] ~ 8"
        self.side_stream_radius = 0.2445

        # reservoir length, injection to production well [m]
        self.res_length = None

        # Number of side streams
        self.n_streams = 1

        # spacing between lateral streams [m]
        self.stream_spacing = 100

        # angle of the lateral wells - horizontal by default [degrees]
        self.angle = 180

        # time [years]
        self.time_years = None

        # surface average ambient air temp [C]
        self.T_surface_air_C = 10

        # surface rock temp [C]
        self.T_surface_rock_C = 10

        # earth's geothermal gradient [C/m], always positive
        self.dT_dz = 0.035

        self.fluid = 'CO2'

        # CO2 properties
        # Specific heat ratio
        self.k = 1.289

        self.pcrit = CP.PropsSI('pcrit', 'CO2')  # Pa
        self.Tcrit = CP.PropsSI('Tcrit', 'CO2') - 273.15  # C

        # configurations [config == 1:4]
        self.config = 1

        # cooling water system
        self.T_cooling_water = 10  # C
        self.P_cooling_water = 5e5  # Pa
        self.eta_cooling_water_pump = 0.75

        # heat exchanger parameters
        self.U_dirty = 600  # W/m^2/K

        # splitter parameters
        self.S_ratio = 1  # Split ratio --> S = m1/m2; mtotal = m1+m2

        # toggle whether ideal S_ratio is calculated for config 3|| default
        # is off.
        self.find_opt_S_ratio_toggle = 0

        self.optimizationMode = 'MinimizeLCOE_Brownfield'
        # optimizationMode = 'MinimizeLCOE_Greenfield';
        # optimizationMode = 'MaximizePower';
        # silicaPrecip = 'PreventSilica';
        self.silicaPrecip = 'IgnoreSilica'
        # coolingMode = 'Wet';
        self.coolingMode = 'Dry'
        self.orcFluid = 'R245fa'
        # orcFluid = 'R600a'
        # thickness for porous heat depletion
        self.thickness = 100
        # wellFieldType
        # (Note: 5spot and 5spot_ManyN with N_5spot=1 are the same thing)
        # 'Doublet','5spot','5spot_SharedNeighbor','5spot_ManyN'
        self.wellFieldType = 'Doublet'
        # N_5spot is only used with the '5spot_ManyN' wellfield type.
        # N_5spot is the sqrt of the number of 5
                
        # O&M Fraction
        self.F_OM = 0.045

        # Discount Rate
        self.discountRate = 0.096

        # Financial Lifetime
        self.Lifetime = 25

        # Capacity Factor
        self.CapacityFactor = 0.9

        # Cost Year
        self.costYear = 2022

        # 2ellCostType. Ideal values are the technological limit, Baseline
        # values are the current state of the art.
        # allowed: 'Ideal','Baseline'
        self.wellCostType = 'Baseline'

        # Success rate of wells
        self.SuccessRate_well = 0.95

        # monitoring well diameter
        self.monitoringWellDiameter = 0.216

        # porous reservoir simulation
        self.modelResPressureTransient = False

        # porous reservoir simulation
        self.modelResTemperatureDepletion = True

        # CPG params
        self.eta_cpg_turbine = 0.78
        self.eta_cpg_pump = 0.9
        self.eta_cpg_compressor = 0.78

        # Turbine parameters
        self.T_turb_out_design = 65  # C
        self.T_turb_out_min = 17

        # Is phase change allowed in the turbine from supercritical
        # to liquid?
        self.PhaseChangeAllowed = False

        # Overall heat transfer coefficient for orc HX
        self.HX_overallHeatTransferCoefficient = 500  # W/m^2-K

        # Pressure drop across a heat exchanger
        self.dP_hex = 0.5e5  # Pa (0.5 bar)

        # ORC params
        self.dT_approach = 7
        self.dT_orc_pinch = 5
        self.eta_orc_pump = 0.9
        self.eta_orc_turbine = 0.8
        self.orcCycleType = 'Subcritical'
        # orcCycleType = 'Supercritical'
        self.orcModel = 'Simulation'
        # orcModel = 'Geophires'

        # hasSurfaceGatheringSystem (true/false). Sets if there is a system of
        # pipes connecting the surface injection and production wells.

        # The injection and production wells are nearby
        self.hasSurfaceGatheringSystem = False

        # water params
        self.PumpDepth = 500
        self.maxPump_dP = 10e6  # 10 MPa max
        self.eta_pump = 0.75

        # well parameters
        self.useWellboreHeatLoss = True
        # These three should be inputs into this function
        self.k_rock = 2.1  # W/m-C
        self.C_rock = 1000  # J/kg-C
        self.rho_rock = 2650  # kg/m^3
        # Friction Factor
        self.epsilon = 55 * 1e-6  # 55 um

        # simulation in steady state 1 = yes, 0 = no (start-up conditions)
        self.steadystate = 1

        # simulatorType 'genGEO' or 'GEOPHIRES'
        self.simulatorType = 'genGEO'

        ## FINANCIAL PARAMETERS
        self.PPI = PPIs.load_PPI(self)