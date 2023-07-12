import numpy as np
from CapitalCost_SurfacePlantCalculator import CapitalCost_SurfacePlantCalculator as SurfacePlantCalc
from CapitalCost_Well import CapitalCost_Well
from CapitalCost_Wellfield import CapitalCost_Wellfield

class CapitalCost_CPG:
    @staticmethod
    def CapitalCost(result_turbine, result_HEX, result_cooling_pump, result_cooling_tower, result_tank, params):
        # Power values are all in kW for a single Injection/Production (IP) pair.

        # Capital Cost is
        # 1. Surface Plant
        # 2. Gathering System
        # 3. Wells
        # 4. Wellfield Development
        # 5. Exploration
        # 6. Well Stimulation

        # 1. Surface Plant
        ######################################################################
        C_TCI_plant = CapitalCost_CPG.SurfacePlant(result_turbine, result_HEX, result_cooling_pump, result_cooling_tower, result_tank, params)

        # 2. Gathering System
        ######################################################################
        if params.hasSurfaceGatheringSystem == True:
            C_gatheringSystem = CapitalCost_CPG.SurfacePipe(params)
        else:
            C_gatheringSystem = 0

        # 3. Well Cost
        ######################################################################
        # Number of wells to pay for

        # If conduction system, the wellLength is the depth plus half reservoir length
        C_well_vertical = CapitalCost_Well(params.depth, 0, np.mean([params.well_radius, params.side_stream_radius]), params)
        C_well_horizontal = CapitalCost_Well(0, params.res_length/2, params.side_stream_radius, params)

        C_well_vertical_green_saved = CapitalCost_Well(3000, 0, params.well_radius, params)
        C_well_vertical_grean_extended = CapitalCost_Well(params.depth, 0, params.well_radius, params) - C_well_vertical_green_saved

        C_wells_horizontal = 2 * C_well_horizontal * params.n_streams

        C_well_brownfield = C_well_vertical + C_well_vertical_grean_extended + C_wells_horizontal
        C_well_greenfield = 2 * C_well_vertical + C_wells_horizontal

        # 4. Wellfield Development
        ######################################################################
        C_wellfield = CapitalCost_Wellfield(params)

        # 5. Exploration
        ######################################################################
        # No exploration costs if conduction system
        C_exploration = 0

        # 6. Well Stimulation
        ######################################################################
        C_stimulation = 0

        C_brownfield = C_TCI_plant + C_gatheringSystem + C_well_brownfield + C_exploration + C_stimulation
        C_greenfield = C_TCI_plant + C_gatheringSystem + C_well_greenfield + C_wellfield + C_exploration + C_stimulation

        # Operational Costs
        # Commodities
        # Liquid CO2 procurement costs
        C_tCO2 = 390  # 2022 USD/MT
        result = {}
        result['CostSurfacePlant']  = C_TCI_plant
        result['C_gatheringSystem'] = C_gatheringSystem
        result['C_wells_production'] =  C_well_vertical
        result['C_wells_horizontal'] = C_wells_horizontal
        result['C_wellfield'] = C_wellfield
        result['C_exploration'] = C_exploration
        result['C_stimulation'] =  C_stimulation
        result['C_brownfield'] =  C_brownfield
        result['C_greenfield'] = C_greenfield

        return result
    
    
    @staticmethod
    def SurfacePlant(result_turbine, result_HEX, result_cooling_pump, result_cooling_tower, result_tank, params):
        # Calculate all the heat exchanger costs
        C_BM_hex = SurfacePlantCalc.HEX(result_HEX, params)

        # Calculate the turbine costs
        if result_turbine['W_turbine'] > 0:
            C_BM_turbine = SurfacePlantCalc.turbine(result_turbine, params)
        else:
            C_BM_turbine = 0

        # Calculate the pumping costs
        if result_cooling_pump['m_dot'] > 0:
            C_BM_cooling_water_pump = SurfacePlantCalc.pump(result_cooling_pump, params)
        else:
            C_BM_cooling_water_pump = 0

        # Calculate the cooling tower costs
        if result_cooling_tower['W_cooling_tower'] > 0:
            C_BM_cooling_tower = SurfacePlantCalc.cooling_tower(result_cooling_tower, params)
        else:
            C_BM_cooling_tower = 0

        # Auxillary and start_up costs
        # Calculate the compressor costs
        # C_BM_compressor = CapitalCost_SurfacePlantfunc.compressor(result_compressor, params)

        # Calculate the storage tanks
        C_BM_storage_tank = SurfacePlantCalc.tank(result_tank, params)

        # Calculate total bare-module costs
        C_TBM = sum(C_BM_hex) + C_BM_turbine + C_BM_cooling_water_pump + C_BM_cooling_tower + C_BM_storage_tank  # USD

        C_TCI = SurfacePlantCalc.TCI(C_TBM)  # USD

        # Operating costs
        # Calculate the Refrigeration costs

        return C_TCI
    
    
