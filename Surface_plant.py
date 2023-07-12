import CoolProp.CoolProp as CP
import numpy as np
from FrictionFactor import FrictionFactor

class Surface_plant:
    @staticmethod
    def plant_dP(P_surface_in, T_surface_in, params):
        # Calculate the enthalpy at the surface inlet temperature
        h_surface_in = CP.PropsSI('HMASS', 'P', P_surface_in, 'T', T_surface_in + 273.15, params.fluid)
        
        # Calculate the density at the surface conditions
        rho_surface = CP.PropsSI('DMASS', 'P', P_surface_in, 'HMASS', h_surface_in, params.fluid)
        
        # Calculate the friction factor
        ff = FrictionFactor(params.well_radius, P_surface_in, h_surface_in, params.m_dot, params)
        
        # Calculate the pressure drop in the surface pipes
        if params.hasSurfaceGatheringSystem == True:
            dP_surfacePipes = ff * params.res_length / (params.well_radius * 2) ** 5 * 8 * params.m_dot ** 2 / rho_surface / np.pi ** 2
        else:
            dP_surfacePipes = 0
        
        # Calculate the surface outlet pressure
        P_surface_out = P_surface_in - dP_surfacePipes
        
        return P_surface_out
    
    @staticmethod
    def TurbWorkfunc(P_turb_in, T_turb_in, params):
        if params.config == 4:
            W_turbine = 0
            T_turb_out_irr = T_turb_in
            P_turb_out_irr = P_turb_in
            s_turb_in = CP.PropsSI('SMASS', 'P', P_turb_in, 'T', T_turb_in + 273.15, params.fluid)
        else:
            h_turb_in = CP.PropsSI('HMASS', 'P', P_turb_in, 'T', T_turb_in + 273.15, params.fluid)
            s_turb_in = CP.PropsSI('SMASS', 'P', P_turb_in, 'T', T_turb_in + 273.15, params.fluid)
            
            if params.config == 2:
                P_turb_out_iter = 115e5
            else:
                P_turb_out_iter = 74e5
            
            toggle = 0
            go = 1
            iter = 1
            
            while go == 1:
                h_turb_out_isen = CP.PropsSI('HMASS', 'P', P_turb_out_iter, 'SMASS', s_turb_in, params.fluid)
                
                h_turb_out_irr = h_turb_in - params.eta_cpg_turbine * (h_turb_in - h_turb_out_isen)
                s_turb_out_irr = CP.PropsSI('SMASS', 'P', P_turb_out_iter, 'HMASS', h_turb_out_irr, params.fluid)
                T_turb_out_irr = CP.PropsSI('T', 'P', P_turb_out_iter, 'SMASS', s_turb_out_irr, params.fluid)
                phase = CP.PhaseSI('P', P_turb_out_iter, 'SMASS', s_turb_in, params.fluid)
                
                P_turb_out_iter = P_turb_out_iter - 0.1e5
                iter = iter + 1
                
                if toggle == 1:
                    go = 0
                elif params.config == 2 or params.config == 5:
                    if T_turb_out_irr < params.T_turb_out_design + 273.15:
                        P_turb_out_iter = P_turb_out_iter + 0.6e5
                        toggle = 1
                else:
                    if phase == 'liquid' or 'twophase' in phase:
                        P_turb_out_iter = P_turb_out_iter + 0.3e5
                        toggle = 1
            
            P_turb_out_irr = P_turb_out_iter
            
            if params.config == 3:
                m_dot_turb = params.m_dot * params.S_ratio
            elif params.config == 4:
                m_dot_turb = 0
            else:
                m_dot_turb = params.m_dot
            
            W_turbine = m_dot_turb * (h_turb_in - h_turb_out_irr)
            
            if W_turbine < 0:
                raise Exception('Turbine Power is Negative')
            
        result = {
        'W_turbine': W_turbine,
        'T': T_turb_out_irr - 273.15,
        'P': P_turb_out_irr,
        's_turb_in': s_turb_in
        }
        return result
    
    @staticmethod
    def Compfunc(P_comp_in, T_comp_in, P_comp_out, T_comp_out, params):
        h_comp_in = CP.PropsSI('HMASS', 'P', P_comp_in, 'T', T_comp_in + 273.15, params.fluid)
        h_comp_out_isen = CP.PropsSI('HMASS', 'P', P_comp_out, 'T', T_comp_out + 273.15, params.fluid)
        h_comp_out_irr = h_comp_in - params.eta_cpg_compressor * (h_comp_in - h_comp_out_isen)

        # The compressor should at start-up conditions pressurize the
        # system until reaching the highest density until the
        # thermosyphon effect is reached; an additional 25% should be
        # added as a design margin.

        W_compressor = h_comp_out_irr - h_comp_in

        result = {'W_compressor': W_compressor}

        return result

    @staticmethod
    def HEXfunc(P_turb_out_irr, T_turb_out_irr, T_well_out, P_well_out, params):
        if params.config == 1:
            P_hex_hot_in = P_turb_out_irr  # Pa
            T_hex_hot_in = T_turb_out_irr  # C
            T_hex_hot_out = CP.PropsSI('T', 'P', P_hex_hot_in - params.dP_hex, 'Q', 0, params.fluid) - 273.15  # J/kg
            T_hex_cold_in = params.T_cooling_water  # C
            T_hex_cold_out = T_hex_hot_in - params.dT_approach  # C

            h_hex_hot_in = CP.PropsSI('HMASS', 'P', P_hex_hot_in, 'T', T_hex_hot_in + 273.15, params.fluid)
            h_hex_hot_satvap = CP.PropsSI('HMASS', 'P', P_hex_hot_in - params.dP_hex, 'Q', 1, params.fluid)  # J/kg
            h_hex_hot_out = CP.PropsSI('HMASS', 'P', P_hex_hot_in - params.dP_hex, 'Q', 0, params.fluid)  # J/kg
            h_hex_cold_in = CP.PropsSI('HMASS', 'P', params.P_cooling_water, 'T', T_hex_cold_in + 273.15, 'water')  # J/kg
            h_hex_cold_out = CP.PropsSI('HMASS', 'P', params.P_cooling_water - params.dP_hex, 'T', T_hex_cold_out + 273.15, 'water')  # J/kg

            m_dot_water_CoolLoop = (params.m_dot * (h_hex_hot_in - h_hex_hot_out)) / (h_hex_cold_out - h_hex_cold_in)

            # proportion of desuperheating from total cooling process.
            cool_frac = (h_hex_hot_in - h_hex_hot_satvap) / (h_hex_hot_in - h_hex_hot_out)

            # number of shells in series
            N_1 = Surface_plant.HEX_shellcount(T_hex_hot_in, T_hex_hot_out, T_hex_cold_in, T_hex_cold_out)
            LMTD_1 = Surface_plant.LMTD(T_hex_hot_in, T_hex_hot_out, T_hex_cold_in, T_hex_cold_out, N_1, cool_frac)  # K
            Q_hex_1 = params.m_dot * (h_hex_hot_in - h_hex_hot_out)  # W
            Q_hex_2 = 0
            Q_net = 0

            result_cooling_tower = Surface_plant.CoolingTowerfunc(params.m_dot, P_hex_hot_in, T_hex_hot_in, h_hex_hot_out, params)

            A_1 = Q_hex_1 / params.U_dirty / LMTD_1  # m^2
            A_2 = 0
            N_2 = 0

            T_CO2_out = T_hex_hot_out
            T_water_out = T_hex_cold_out
            P_CO2_out = P_hex_hot_in - params.dP_hex
            P_water_out = params.P_cooling_water - params.dP_hex
            m_dot_water = 0
            P_tubeside_1 = P_hex_hot_in - params.dP_hex / 2
            P_tubeside_2 = 0

        elif params['config'] == 2 or params['config'] == 4:
            if params['config'] == 2:
                T_hex_hot_in = T_turb_out_irr  # C
                P_hex_hot_in = P_turb_out_irr  # Pa
            elif params['config'] == 4:
                T_hex_hot_in = T_well_out  # C
                P_hex_hot_in = P_well_out  # Pa
            
            T_hex_hot_out = params['T_cooling_water'] + params['dT_approach']  # C
            T_hex_cold_in = params['T_cooling_water']  # C
            T_hex_cold_out = 50  # C
            
            # number of shells in series
            N_1 = Surface_plant.HEX_shellcount(T_hex_hot_in, T_hex_hot_out, T_hex_cold_in, T_hex_cold_out)
            
            LMTD_1 = Surface_plant.LMTD(T_hex_hot_in, T_hex_hot_out, T_hex_cold_in, T_hex_cold_out, N_1, 0)  # K
            
            h_hex_hot_in = CP.PropsSI('HMASS', 'P', P_hex_hot_in, 'T', T_hex_hot_in + 273.15, params['fluid'])
            h_hex_hot_out = CP.PropsSI('HMASS', 'P', P_hex_hot_in - params['dP_hex'], 'T', T_hex_hot_out + 273.15 - 0.001, params['fluid'])  # J/kg
            h_hex_cold_in = CP.PropsSI('HMASS', 'P', params['P_cooling_water'], 'T', T_hex_cold_in + 273.15, 'water')  # J/kg
            h_hex_cold_out = CP.PropsSI('HMASS', 'P', params['P_cooling_water'] - params['dP_hex'], 'T', T_hex_cold_out + 273.15, 'water')  # J/kg
            
            Q_hex_1 = params['m_dot'] * (h_hex_hot_in - h_hex_hot_out)  # W
            Q_hex_2 = 0
            Q_net = Q_hex_1
            
            result_cooling_tower = {}

            result_cooling_tower['Q_desuperheating'] = 0  # W
            result_cooling_tower['Q_condensing'] = 0  # W
            result_cooling_tower['Q_cooling_tower'] = 0  # W
            result_cooling_tower['W_cooling_tower'] = 0  # W
            result_cooling_tower['dT_range'] = 0  # K
            
            A_1 = Q_hex_1 / params['U_dirty'] / LMTD_1  # m^2
            A_2 = 0
            N_2 = 0
            
            m_dot_water = Q_hex_1 / (h_hex_cold_out - h_hex_cold_in)
            m_dot_water_CoolLoop = 0
            T_CO2_out = T_hex_hot_out
            T_water_out = T_hex_cold_out
            P_CO2_out = P_hex_hot_in - params['dP_hex']
            P_water_out = params['P_cooling_water'] - params['dP_hex']
            P_tubeside_1 = P_hex_hot_in - params['dP_hex'] / 2
            P_tubeside_2 = 0
        elif params.config == 3:
            m_dot_sco2_1 = params.m_dot * params.S_ratio
            m_dot_sco2_2 = params.m_dot * (1 - params.S_ratio)

            # HEX stream 1
            P_hex1_hot_in = P_turb_out_irr  # Pa
            T_hex1_hot_in = T_turb_out_irr  # C
            T_hex1_cold_in = params.T_cooling_water  # C
            T_hex1_cold_out = T_hex1_hot_in - params.dT_approach

            h_hex1_hot_in = CP.PropsSI('HMASS', 'P', P_hex1_hot_in, 'T', T_hex1_hot_in + 273.15, params.fluid)  # J/kg
            h_hex1_cold_in = CP.PropsSI('HMASS', 'P', params.P_cooling_water, 'T', T_hex1_cold_in + 273.15, 'water')  # J/kg
            h_hex1_cold_out = CP.PropsSI('HMASS', 'P', params.P_cooling_water - params.dP_hex, 'T', T_hex1_cold_out + 273.15, 'water')  # J/kg

            # HEX stream 2
            P_hex2_hot_in = P_well_out  # Pa
            T_hex2_cold_in = T_hex1_cold_out
            T_hex2_hot_in = T_well_out  # C
            T_hex2_hot_out = T_hex2_cold_in + params.dT_approach  # C
            T_hex2_cold_out = 50  # C

            h_hex2_hot_in = CP.PropsSI('HMASS', 'P', P_hex2_hot_in, 'T', T_hex2_hot_in + 273.15, params.fluid)  # J/kg
            h_hex2_hot_out = CP.PropsSI('HMASS', 'P', P_hex2_hot_in - params.dP_hex, 'T', T_hex2_hot_out + 273.15, params.fluid)  # J/kg
            h_hex2_cold_out = CP.PropsSI('HMASS', 'P', params.P_cooling_water - 2 * params.dP_hex, 'T', T_hex2_cold_out + 273.15, 'water')  # J/kg

            # Combined streams
            h_hex1_hot_satvap = CP.PropsSI('HMASS', 'P', P_hex1_hot_in, 'Q', 1, params.fluid)  # J/kg
            h_hex_hot_satliq = CP.PropsSI('HMASS', 'P', P_hex1_hot_in - params.dP_hex, 'Q', 0, params.fluid)  # J/kg

            P_CO2_out = P_hex1_hot_in - params.dP_hex
            T_hot_satliq = CP.PropsSI('T', 'P', P_CO2_out, 'HMASS', h_hex_hot_satliq, params.fluid) - 273.15  # C
            
            if params.find_opt_S_ratio_toggle == 0:
                h_hex1_hot_out = (m_dot_sco2_2 / m_dot_sco2_1) * (h_hex_hot_satliq - h_hex2_hot_out) + h_hex_hot_satliq
                if m_dot_sco2_1 > 0:
                    T_hex1_hot_out = CP.PropsSI('T', 'P', P_CO2_out, 'HMASS', h_hex1_hot_out, params.fluid) - 273.15  # C
                T_CO2_out = T_hot_satliq
                
                if h_hex1_hot_out > h_hex1_hot_in:
                    cooling_required = 0
                else:
                    cooling_required = 1

            elif params.find_opt_S_ratio_toggle == 1:
                cooling_required = 0
                
                params.S_ratio = (h_hex_hot_satliq - h_hex2_hot_out) / (h_hex1_hot_in - h_hex2_hot_out)
                m_dot_sco2_1 = params.m_dot * params.S_ratio
                m_dot_sco2_2 = params.m_dot * (1 - params.S_ratio)
                
                T_hex2_cold_in = T_hex1_cold_in
                h_hex2_cold_in = h_hex1_cold_in
                m_dot_water = m_dot_sco2_2 * (h_hex2_hot_in - h_hex2_hot_out) / (h_hex2_cold_out - h_hex2_cold_in)
                m_dot_water_CoolLoop = 0

            if cooling_required == 0:
                h_hex_hot_combined_out = params.S_ratio * h_hex1_hot_in + (1 - params.S_ratio) * h_hex2_hot_out
                T_CO2_out = CP.PropsSI('T', 'P', P_CO2_out, 'HMASS', h_hex_hot_combined_out, params.fluid) - 273.15  # C

            if params.find_opt_S_ratio_toggle == 1:
                pass

            elif cooling_required == 0:
                h_hex2_cold_in = h_hex1_cold_in
                T_hex2_cold_in = T_hex1_cold_in
                m_dot_water = (m_dot_sco2_2 * (h_hex2_hot_in - h_hex2_hot_out)) / (h_hex2_cold_out - h_hex2_cold_in)
                m_dot_water_CoolLoop = 0

            elif cooling_required == 1:
                m_dot_water_hex_1 = m_dot_sco2_1 * (h_hex1_hot_in - h_hex1_hot_out) / (h_hex1_cold_out - h_hex1_cold_in)
                h_hex2_cold_in = h_hex1_cold_out
                T_hex2_cold_in = CP.PropsSI('T', 'P', params.P_cooling_water - params.dP_hex, 'HMASS', h_hex2_cold_in, 'water') - 273.15  # C
                m_dot_water = (m_dot_sco2_2 * (h_hex2_hot_in - h_hex2_hot_out)) / (h_hex2_cold_out - h_hex2_cold_in)
                m_dot_water_CoolLoop = m_dot_water_hex_1 - m_dot_water

                cool_frac_1 = (h_hex1_hot_in - h_hex1_hot_satvap) / (h_hex1_hot_in - h_hex1_hot_out)

            if params.find_opt_S_ratio_toggle == 0 and cooling_required == 1:
                N_1 = Surface_plant.HEX_shellcount(T_hex1_hot_in, T_hex1_hot_out, T_hex1_cold_in, T_hex1_cold_out)
                LMTD_1 = Surface_plant.LMTD(T_hex1_hot_in, T_hex1_hot_out, T_hex1_cold_in, T_hex1_cold_out, N_1, cool_frac_1)  # K
                Q_hex_1 = m_dot_sco2_1 * (h_hex1_hot_in - h_hex1_hot_out)  # W
                A_1 = Q_hex_1 / params.U_dirty / LMTD_1  # m^2
                P_water_out = params.P_cooling_water - params.dP_hex
                P_tubeside_1 = P_hex1_hot_in - params.dP_hex / 2

                result_cooling_tower = Surface_plant.CoolingTowerfunc(m_dot_sco2_1, P_hex1_hot_in, T_hex1_hot_in, h_hex1_hot_out, params)

            else:
                N_1 = 0
                Q_hex_1 = 0  # W

                result_cooling_tower = {}
                result_cooling_tower['Q_desuperheating'] = 0  # W
                result_cooling_tower['Q_condensing'] = 0  # W
                result_cooling_tower['Q_cooling_tower'] = 0  # W
                result_cooling_tower['W_cooling_tower'] = 0  # W
                result_cooling_tower['dT_range'] = 0  # K

                A_1 = 0  # m^2
                P_water_out = params.P_cooling_water
                P_tubeside_1 = 0

            # Heat Exchanger 2 (large)
            # number of shells in series
            N_2 = Surface_plant.HEX_shellcount(T_hex2_hot_in, T_hex2_hot_out, T_hex2_cold_in, T_hex2_cold_out)

            LMTD_2 = Surface_plant.LMTD(T_hex2_hot_in, T_hex2_hot_out, T_hex2_cold_in, T_hex2_cold_out, N_2, 0)  # K
            Q_hex_2 = m_dot_sco2_2 * (h_hex2_hot_in - h_hex2_hot_out)  # W
            Q_net = Q_hex_1 + Q_hex_2
            A_2 = Q_hex_2 / params.U_dirty / LMTD_2  # m^2

            T_water_out = T_hex2_cold_out
            P_CO2_out = P_hex1_hot_in - params.dP_hex
            P_water_out = P_water_out - params.dP_hex
            P_tubeside_2 = P_hex2_hot_in - params.dP_hex / 2
            
        result = {}

        result['Q_hex_1'] = Q_hex_1
        result['Q_hex_2'] = Q_hex_2
        result['result_cooling_tower'] = result_cooling_tower

        # increasing the surface area by 50% to account for fluctuations in cooling water inlet temperature
        # and thus reduced heat exchange efficiency
        result['Q_net'] = Q_net
        result['A_1'] = A_1 * 1.5
        result['A_2'] = A_2 * 1.5
        result['N_1'] = N_1
        result['N_2'] = N_2
        result['m_dot_water'] = m_dot_water
        result['m_dot_water_CoolLoop'] = m_dot_water_CoolLoop
        result['T_water_out'] = T_water_out
        result['T_CO2_out'] = T_CO2_out
        result['P_water_out'] = P_water_out
        result['P_CO2_out'] = P_CO2_out
        result['P_tubeside_1'] = P_tubeside_1
        result['P_tubeside_2'] = P_tubeside_2
        return result

    @staticmethod
    def LMTD(T_hot_in, T_hot_out, T_cold_in, T_cold_out, N, cool_frac):
        dT_in = abs(T_hot_in - T_cold_out)  # K
        dT_out = abs(T_hot_out - T_cold_in)  # K
        LMTD = (dT_in - dT_out) / np.log(dT_in / dT_out)

        R = (T_hot_in - T_hot_out) / (T_cold_out - T_cold_in)
        P = (T_cold_out - T_cold_in) / (T_hot_in - T_cold_in)

        S = (R**2 + 1)**0.5 / (R - 1)
        W = ((1 - P * R) / (1 - P))**(1 / N)
        F_T = (S * np.log(W)) / (np.log((1 + W - S + S * W) / (1 + W + S - S * W)))

        if T_hot_in == T_hot_out:
            pass
        elif cool_frac > 0:
            LMTD = LMTD * (cool_frac * F_T + (1 - cool_frac))
        else:
            LMTD = LMTD * F_T

        return LMTD

    @staticmethod
    def HEX_shellcount(T_hot_in, T_hot_out, T_cold_in, T_cold_out):
        m_cold = T_cold_out - T_cold_in
        m_hot = T_hot_in - T_hot_out

        y = T_hot_out
        N = 0
        x = 0
        while x < 0.95:
            x = (y - T_cold_in) / m_cold
            y = m_hot * x + T_hot_out
            N = N + 1

        return N
    
    @staticmethod
    def ParasiticPowerFraction_CoolingTower(T_ambient_C, dT_approach_CT, dT_range_CT, coolingMode):
        if coolingMode == 'Wet':
            # wet
            a_cool_wet = 1.20e0
            b_cool_wet = 0
            c_cool_wet = -3.79e-3
            d_cool_wet = 1.95e-2
            f_cooling = a_cool_wet * (1 / dT_approach_CT) + b_cool_wet * (T_ambient_C + 273.15) + c_cool_wet * (T_ambient_C + 273.15) / dT_approach_CT + d_cool_wet * (1 / (dT_approach_CT + dT_range_CT))
            
            a_cond_wet = 1.65e0
            b_cond_wet = -6.24e-6
            c_cond_wet = -5.03e-3
            d_cond_wet = 0
            f_condensing = a_cond_wet * (1 / dT_approach_CT) + b_cond_wet * (T_ambient_C + 273.15) + c_cond_wet * (T_ambient_C + 273.15) / dT_approach_CT + d_cond_wet * (1 / (dT_approach_CT + dT_range_CT))
        
        elif coolingMode == 'Dry':
            # dry
            a_cool_dry = 7.65e-1
            b_cool_dry = 0
            c_cool_dry = 0
            d_cool_dry = 1.28e-1
            f_cooling = a_cool_dry * (1 / dT_approach_CT) + b_cool_dry * (T_ambient_C + 273.15) + c_cool_dry * (T_ambient_C + 273.15) / dT_approach_CT + d_cool_dry * (1 / (dT_approach_CT + dT_range_CT))
        
            a_cond_dry = 6.19e-1
            b_cond_dry = 0
            c_cond_dry = 0
            d_cond_dry = 0
            f_condensing = a_cond_dry * (1 / dT_approach_CT) + b_cond_dry * (T_ambient_C + 273.15) + c_cond_dry * (T_ambient_C + 273.15) / dT_approach_CT + d_cond_dry * (1 / (dT_approach_CT + dT_range_CT))
        else:
            raise Exception('Unknown Cooling Mode')
        
        return f_cooling, f_condensing

    def CoolingTowerfunc(m_dot, P_cond_in, T_cond_in, h_cond_out, params):
        # heat rejection
        h_cond_in = CP.PropsSI('HMASS', 'P', P_cond_in, 'T', T_cond_in + 273.15, params.fluid)

        if P_cond_in < params.pcrit:
            h_satVapor = CP.PropsSI('HMASS', 'P', P_cond_in, 'Q', 1, params.fluid)
        else:
            h_satVapor = 1e7

        T_cond_out = CP.PropsSI('T', 'P', P_cond_in - params.dP_hex, 'HMASS', h_cond_out, params.fluid) - 273.15

        if h_cond_out > h_satVapor or P_cond_in > params.pcrit:
            # only desuperheating needed, no condensation required
            Q_cooler_part = m_dot * (h_cond_in - h_cond_out)
            Q_condenser_part = 0
            dT_range = T_cond_in - T_cond_out
        elif h_cond_in > h_satVapor:
            # desuperheating needed
            Q_cooler_part = m_dot * (h_cond_in - h_satVapor)
            Q_condenser_part = m_dot * (h_satVapor - h_cond_out)
            dT_range = T_cond_in - T_cond_out
        else:
            # no desuperheating
            Q_cooler_part = 0
            Q_condenser_part = m_dot * (h_cond_in - h_cond_out)
            dT_range = 0
        
        f_cooling, f_condensing = Surface_plant.ParasiticPowerFraction_CoolingTower(params.T_surface_air_C, params.dT_approach, dT_range, params.coolingMode)
        W_cooler_part = f_cooling * Q_cooler_part
        W_condenser_part = f_condensing * Q_condenser_part
        Q_cooling_tower = Q_cooler_part + Q_condenser_part 
        W_cooling_tower = W_cooler_part + W_condenser_part
        
        result = {
            'Q_desuperheating': Q_cooler_part,
            'Q_condensing': Q_condenser_part,
            'Q_cooling_tower': Q_cooling_tower,
            'W_cooling_tower': W_cooling_tower,
            'dT_range': dT_range
        }
        
        return result
    
    @staticmethod
    def PumpWorkfunc(dP, m_dot, fluid, params):
        # Calculate net pump power value
        rho = CP.PropsSI('DMASS', 'P', params.P_cooling_water, 'T', params.T_cooling_water + 273.15, fluid)  # kg/m^3
        V_dot = m_dot / rho

        W_pump = dP * V_dot / params.eta_cooling_water_pump

        result = {'W_pump': W_pump, 'dP': dP, 'm_dot': m_dot, 'rho': rho}

        return result
    
    @staticmethod
    def tank(injWell, reservoir, prodWell, turbine, params):

        m1 = np.pi * (np.tile(injWell['wellRadius'], (params.n_streams, 1))).T ** 2 * injWell['Density'] * injWell['dL']
        m2 = np.pi * (np.tile(reservoir['wellRadius'], (params.n_streams, 1))).T ** 2 * reservoir['Density'] * reservoir['dL'] * np.ones((1, params.n_streams))
        m3 = np.pi * (np.tile(prodWell['wellRadius'], (params.n_streams, 1))).T ** 2 * prodWell['Density']* prodWell['dL']

        m_total = np.sum(m1[:, -1]) + np.sum(m2) + np.sum(m3[:, -1])  # kg

        # store as liquid at P > P(critical)
        # P_store = CP.PropsSI('CO2', 'pcrit') + 0.1e5  # Pa
        P_design_margin = 1.8e5
        P_store = turbine['P'] + P_design_margin  # Pa
        T_design_margin = 5  # C
        T_store = params.T_surface_air_C + T_design_margin  # C

        # Rho at storage condition: P above critical to avoid phase change and rapid decompression
        rho_store = CP.PropsSI('DMASS', 'P', P_store, 'T', T_store + 273.15, params.fluid)  # kg/m^3

        # Store 10% for contingency
        V_total_surface = m_total / rho_store  # m^3
        V_store = V_total_surface * 0.1

        result = {'P_store': P_store, 'm_total': m_total, 'V_total_surface': V_total_surface, 'V_store': V_store}

        return result
    
    @staticmethod
    def NetWorkfunc(W_turbine, W_cooling_pump, W_cooling_tower, W_CO2_SU_comp, params):
        if params.steadystate == 1:
            # calculate net power values
            W_net = W_turbine - W_cooling_pump - W_cooling_tower
        else:
            W_net = W_turbine - W_cooling_pump - W_cooling_tower - W_CO2_SU_comp

        return W_net




