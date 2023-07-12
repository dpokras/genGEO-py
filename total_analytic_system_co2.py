import numpy as np
import CoolProp.CoolProp as CP
from semi_analytic_well import semi_analytic_well
from Surface_plant import Surface_plant
from CapitalCost_CPG import CapitalCost_CPG as CC
from LCOE_Simple import LCOE_Simple

def total_analytic_system_co2(params):
    if params.fluid != 'CO2':
        raise Exception('total_analytic_system_CO2:Wrong_Fluid - Wrong Fluid -- Not CO2')

    T_inj_surface = 25  # C
    
    P_inj_surface = 65e5 # Pa

    P_reservoir = params.depth * 1000 * 9.81
    P_reservoir_max = params.depth * 2500 * 9.81

    dP_threshold = 2e4  # 20 kPa
    iter = 1

    while True:

        # enthalpy state 1 (injection)
        H1 = np.array([CP.PropsSI('HMASS', 'T', T_inj_surface+273.15, 'P', P_inj_surface, params.fluid) * params.m_dot] * params.n_streams)
        h1 = H1 / np.array(params.m_dot)
        P1 = np.array([P_inj_surface] * params.n_streams)

        depths = params.depth - params.thickness * np.linspace(params.n_streams - 1, 0, params.n_streams)
        temp_at_depth = params.T_surface_rock_C + params.dT_dz * depths  # C
        
        result_injWell = semi_analytic_well(P_inj_surface, T_inj_surface, params.T_surface_rock_C, -1*depths, 0, params.m_dot, params)
        result_reservoir = semi_analytic_well(result_injWell['EndPressure'], result_injWell['EndTemp'], temp_at_depth, 0, params.res_length, params.m_dot/params.n_streams, params)

        if result_reservoir['EndPressure'][-1] >= P_reservoir_max:
            raise Exception('TotalAnalyticSystemCO2:ExceedsMaxReservoirPressure', f'Exceeds Max Reservoir Pressure of {P_reservoir_max/1e6:.3f} MPa!')

        result_prodWell = semi_analytic_well(result_reservoir['EndPressure'], result_reservoir['EndTemp'], temp_at_depth, depths, 0, params.m_dot, params)
        
        P_prod_surface_flashed = np.array([result_prodWell['EndPressure'][0]] * params.n_streams)
        T_prod_surface_flashed = CP.PropsSI('T', 'P', P_prod_surface_flashed, 'HMASS', result_prodWell['EndEnthalpy'], params.fluid) - 273.15
        T_prod_surface = np.mean(T_prod_surface_flashed)
        P_prod_surface = np.mean(P_prod_surface_flashed)

        # enthalpy state 2 (injection well)
        H2 = result_injWell['Enthalpy'] * params.m_dot
        h2 = result_injWell['Enthalpy']
        P2 = result_injWell['Pressure']

        # enthalpy state 3 (reservoir)
        H3 = result_reservoir['Enthalpy'] * params.m_dot
        h3 = result_reservoir['Enthalpy']
        P3 = result_reservoir['Pressure']

        # enthalpy state 4 (production well)
        H4 = result_prodWell['Enthalpy'] * params.m_dot
        h4 = result_prodWell['Enthalpy']
        P4 = result_prodWell['Pressure']
        H4_F = np.array([CP.PropsSI('HMASS', 'T', T_prod_surface+273.15, 'P', P_prod_surface, params.fluid) * params.m_dot] * params.n_streams)
        h4_F = H4_F / np.array(params.m_dot)
        P4_F = np.array([P_prod_surface] * params.n_streams)

        # enthalpy state 5 (stream split; S1/2)
        H5S1 = np.array([CP.PropsSI('HMASS', 'T', T_prod_surface+273.15, 'P', P_prod_surface, params.fluid) * params.m_dot * params.S_ratio] * params.n_streams)
        h5S1 = H5S1 / (np.array(params.m_dot) * params.S_ratio)
        H5S2 = np.array([CP.PropsSI('HMASS', 'T', T_prod_surface+273.15, 'P', P_prod_surface, params.fluid) * params.m_dot * (1 - params.S_ratio)] * params.n_streams)
        h5S2 = H5S2 / (np.array(params.m_dot) * (1 - params.S_ratio))
        P5S1 = np.array([P_prod_surface] * params.n_streams)
        P5S2 = np.array([P_prod_surface] * params.n_streams)

        P_surfacePlant = Surface_plant.plant_dP(P_prod_surface, T_prod_surface, params)

        result_turbine = Surface_plant.TurbWorkfunc(P_surfacePlant, T_prod_surface, params)
        W_turbine = result_turbine['W_turbine']

        # enthalpy state 6 (turbine outlet)
        H6 = np.array([CP.PropsSI('HMASS', 'T', result_turbine['T']+273.15, 'P', result_turbine['P'], params.fluid) * params.m_dot * params.S_ratio] * params.n_streams)
        h6 = H6 / (np.array(params.m_dot) * params.S_ratio)
        P6 = np.array([result_turbine['P']] * params.n_streams)

        result_HEX = Surface_plant.HEXfunc(result_turbine['P'], result_turbine['T'], T_prod_surface, P_prod_surface, params)
        Q_net = result_HEX['Q_net']

        # enthalpy state 7 (S1 HEX1 outlet)
        H7 = H6 - np.array(result_HEX['Q_hex_1'])
        h7 = H7 / (params.m_dot * params.S_ratio)
        P7 = np.array([result_HEX['P_tubeside_1']] * params.n_streams)

        H8 = H5S2 - np.array(result_HEX['Q_hex_2'])
        h8 = H8 / (np.array(params.m_dot) * (1 - params.S_ratio))
        P8 = np.array([result_HEX['P_tubeside_2']] * params.n_streams)

        H9 = H8  # pressure is different P9 = P7
        h9 = h8
        P9 = P7

        H10 = np.array([CP.PropsSI('HMASS', 'T', result_HEX['T_CO2_out']+273.15-0.001, 'P', result_HEX['P_CO2_out'], params.fluid) * params.m_dot] * params.n_streams)
        h10 = H10 / np.array(params.m_dot)
        P10 = P7

        # Concat all Hs and Ps
        HS1 = np.concatenate((H1, H2, H3, H4, H4_F, H5S1, H6, H7, H10),1)
        hS1 = np.concatenate((h1, h2, h3, h4, h4_F, h5S1, h6, h7, h10),1)
        PS1 = np.concatenate((P1, P2, P3, P4, P4_F, P5S1, P6, P7, P10),1)

        HS2 = np.concatenate((H1, H2, H3, H4, H4_F, H5S2, H8, H9, H10),1)
        hS2 = np.concatenate((h1, h2, h3, h4, h4_F, h5S2, h8, h9, h10),1)
        PS2 = np.concatenate((P1, P2, P3, P4, P4_F, P5S2, P8, P9, P10),1)

        dP_cooling_water = 2e5  # Pa
        result_cooling_pump = Surface_plant.PumpWorkfunc(dP_cooling_water, result_HEX['m_dot_water_CoolLoop'], 'water', params)
        W_cooling_pump = result_cooling_pump['W_pump']

        W_compressor = 0

        result_cooling_tower = result_HEX['result_cooling_tower']
        W_cooling_tower = result_cooling_tower['W_cooling_tower']
        
        W_net = Surface_plant.NetWorkfunc(W_turbine, W_cooling_pump, W_cooling_tower, W_compressor, params)

        result_tank = Surface_plant.tank(result_injWell, result_reservoir, result_prodWell, result_turbine, params)

        dP_inj = result_turbine['P'] - P_inj_surface
        if dP_inj < dP_threshold:
            break

        P_inj_surface = result_turbine['P']
        T_inj_surface = result_HEX['T_CO2_out']
        iter += 1

    Q_fluid = (result_injWell['Heat'] + result_reservoir['Heat'] + result_prodWell['Heat'])

    result_capitalCost = CC.CapitalCost(result_turbine, result_HEX, result_cooling_pump, result_cooling_tower, result_tank, params)

    result_brownfield = LCOE_Simple(result_capitalCost['C_brownfield'], W_net, Q_net, params)
    result_greenfield = LCOE_Simple(result_capitalCost['C_greenfield'], W_net, Q_net, params)

    result = {}
    result['injection'] = result_injWell
    result['reservoir'] = result_reservoir
    result['productionLower'] = []
    result['productionUpper'] = result_prodWell

    result['result_turbine'] = result_turbine
    result['result_HEX'] = result_HEX
    result['result_cooling_pump'] = result_cooling_pump
    result['result_cooling_tower'] = result_cooling_tower
    result['W_turbine'] = W_turbine
    result['W_net'] = W_net
    result['Q_net'] = Q_net
    result['s_turb_in'] = result_turbine['s_turb_in']

    result['temp_at_depth'] = temp_at_depth
    result['P_reservoir'] = P_reservoir
    result['P_reservoir_max'] = P_reservoir_max
    result['T_prod_surface'] = T_prod_surface
    result['P_prod_surface'] = P_prod_surface
    result['Q_fluid'] = Q_fluid
    result['max_speed'] = params.m_dot / (np.pi * (params.side_stream_radius) ** 2) / result_prodWell['EndDensity'][0]

    result['CapitalCost'] = result_capitalCost
    result['SpCC_W_brownfield'] = result_brownfield['SpCC_W']
    result['SpCC_W_greenfield'] = result_greenfield['SpCC_W']
    result['SpCC_Q_brownfield'] = result_brownfield['SpCC_Q']
    result['SpCC_Q_greenfield'] = result_greenfield['SpCC_Q']
    result['LCOE_brownfield'] = result_brownfield['LCOE']
    result['LCOE_greenfield'] = result_greenfield['LCOE']

    # PH diagram
    result['HS1'] = HS1
    result['HS2'] = HS2
    result['hS1'] = hS1
    result['hS2'] = hS2
    result['PS1'] = PS1
    result['PS2'] = PS2
    return result
