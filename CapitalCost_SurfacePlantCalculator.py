import numpy as np
import CoolProp as CP
from PPIs import PPIs

class CapitalCost_SurfacePlantCalculator:
    @staticmethod
    def HEX(input, params):
        # STEADY STATE COSTS

        # USD to EURO Conversion (May 2023)
        # 1 USD = 0.92 EURO

        # Heat exchangers
        # Description: Floating head type Shell & tube HEX. Carbon steel except
        # tube-side Cr-Mo steel for corrosion resistance. sCO2 is tube-side;
        # H2O is shell-side.

        # REFERENCE: Product and Process Design Principles: Synthesis,
        # Analysis, and Evaluation 4th Ed.

        # Heat exchanger 1

        # reference CEPCI; 2013 (average)
        PPI_ref = PPIs.return_PPI('PPI_HX', 2013, params)  # 2013 330.167
        PPI_current = PPIs.return_PPI('PPI_HX', params.costYear, params)  # 2022 494.764

        # inputs
        area = [input['A_1'], input['A_2']]
        area = [x * 10.764 for x in area]  # ft^2
        P = [input['P_tubeside_1'], input['P_tubeside_2']]

        # free-of-board costs
        C_B = [np.exp(12.0310 - 0.8709 * np.log(x) + 0.09005 * (np.log(x)) ** 2) if x != 0 else 0 for x in area]  # USD

        # material factor -- S/T-side = carbon steel/Cr-Mo steel
        a = 1.55
        b = 0.05
        F_M = [a + (x / 100) ** b for x in area]

        # length factor L = 6069mm
        F_L = 1.0

        # tubeside pressure factor - cost correlation (P > 10.342 bar)
        # REFERENCE: CAMARAZA-MEDINA et al. (2021) https://doi.org/10.52292/j.laar.2021.713
        F_P = []
        ID = 0.65  # m
        for p in P:
            if p > 10.342e5:  # Pa
                factor = 1 + ((p / 1e5) / 10.342 - 1) * (0.0035 - 0.022 * (ID - 0.3048))  # P in Pa, ID in m
            else:
                factor = 1.0
            F_P.append(factor)

        C_P = [f_p * f_m * F_L * c_b for f_p, f_m, c_b in zip(F_P, F_M, C_B)]  # USD
        C_P = [x * (PPI_current / PPI_ref) for x in C_P]  # USD

        # Guthrie Bare-Module Factor for STHEs
        F_BM = 3.17
        C_BM = [x * F_BM for x in C_P]

        return C_BM
    
    @staticmethod
    def pump(input, params):
        # Pump Cost Estimation
        
        # reference CEPCI; 2013 (average)
        PPI_ref = PPIs.return_PPI('PPI_Pump&Comp', 2013, params)  # 2013 140.167
        PPI_current = PPIs.return_PPI('PPI_Pump&Comp', params.costYear, params)  # 2022 192.573
        
        Q = input['m_dot'] / input['rho'] * 15850  # gallons per minute
        H = (input['dP'] / 1e5) * 10.197 * 3.281  # ft
        
        S = Q * (H) ** 0.5  # Size factor
        
        C_B_pump = np.exp(12.1656 - 1.1448 * np.log(S) + 0.0862 * (np.log(S)) ** 2)  # USD
        
        # pump type factor; 1 stage, 1800rpm, 250 - 5000gpm range, 50 - 500ft head range, 250Hp max motor power
        F_T_pump = 2.0
        
        # material factor; cast iron
        F_M_pump = 1.0
        
        # Electric motor calculations for pump
        mu_P = -0.316 + 0.24015 * (np.log(Q)) - 0.01199 * (np.log(Q)) ** 2
        P_B = Q * H * (input['rho'] / 119.8) / 33000 / mu_P  # Hp
        mu_M = 0.80 + 0.0319 * (np.log(P_B)) - 0.00182 * (np.log(P_B)) ** 2
        P_c = P_B / mu_M  # Hp
        
        C_B_driver = np.exp(5.9332 + 0.16829 * np.log(P_c) - 0.110056 * (np.log(P_c)) ** 2 +
                            0.071413 * (np.log(P_c)) ** 3 - 0.0063788 * (np.log(P_c)) ** 4)  # USD
        
        F_T_driver = 1.3  # 1800rpm enclosed, fan-cooled, 1 to  250 Hp
        
        C_P = C_B_pump * F_T_pump * F_M_pump + C_B_driver * F_T_driver  # USD
        C_P = C_P * (PPI_current / PPI_ref)  # USD
        
        # Guthrie Bare-Module Factor for pumps
        F_BM = 3.30
        C_BM = C_P * F_BM
        
        return C_BM
    
    @staticmethod
    def compressor(input, params):
        PPI_ref = PPIs.return_PPI('PPI_Pump&Comp', 2013, params)  # 2013 140.167
        PPI_current = PPIs.return_PPI('PPI_Pump&Comp', params.costYear, params)  # 2022 192.573

        P_c = input['P_c'] / 745.7  # Hp
        C_B = np.exp(9.1553 + 0.63 * np.log(P_c))  # USD

        F_D = 1.0
        F_M = 2.5

        C_P = C_B * F_M * F_D  # USD
        C_P = C_P * (PPI_current / PPI_ref)  # USD

        F_BM = 2.15  # Guthrie Bare-Module Factor for gas compressors
        C_BM = C_P * F_BM

        return C_BM
    @staticmethod
    def turbine(input, params):
        PPI_ref = PPIs.return_PPI('PPI_T-G', 2003, params)  # 2003 154.033
        PPI_current = PPIs.return_PPI('PPI_T-G', params.costYear, params)  # 2022 245.235

        S_T_fluid = 1.20  # Regular fluid (CO2)
        C_P = 0.67 * (S_T_fluid * 2830 * (input['W_turbine'] / 1e3) ** 0.745 + 3680 * (input['W_turbine'] / 1e3) ** 0.617)

        C_P = C_P * (PPI_current / PPI_ref)  # USD

        F_BM = 1.5  # Guthrie Bare-Module Factor for gas-driven turbines
        C_BM = C_P * F_BM

        return C_BM
    @staticmethod
    def cooling_tower(input, params):
        PPI_ref = PPIs.return_PPI('PPI_Cool', 2019, params)  # 2013 131.650
        PPI_current = PPIs.return_PPI('PPI_Cool', params.costYear, params)  # 2022 272.893

        TDC = 1.2  # Tower Design Coefficient - CO2: 1.2

        if params.coolingMode == 'Wet':
            a_cool = 5.58e3
            b_cool = 0
            c_cool = -1.77e1
            d_cool = 1.96e2
            c_cooling = a_cool * (1 / params.dT_approach) + b_cool * (params.T_surface_air_C + 273.15) + c_cool * (params.T_surface_air_C + 273.15) / params.dT_approach + d_cool * (1 / (params.dT_approach + input['dT_range']))
            a_cond = 4.08e3
            b_cond = -1.54e-2
            c_cond = -1.24e1
            d_cond = 0
            c_condensing = a_cond * (1 / params.dT_approach) + b_cond * (params.T_surface_air_C + 273.15) + c_cond * (params.T_surface_air_C + 273.15) / params.dT_approach + d_cond * (1 / (params.dT_approach + input['dT_range']))
        elif params.coolingMode == 'Dry':
            a_cool = 7.31e3
            b_cool = 0
            c_cool = 0
            d_cool = 1.23e3
            c_cooling = a_cool * (1 / params.dT_approach) + b_cool * (params.T_surface_air_C + 273.15) + c_cool * (params.T_surface_air_C + 273.15) / params.dT_approach + d_cool * (1 / (params.dT_approach + input['dT_range']))
            a_cond = 1.91e3
            b_cond = 0
            c_cond = 0
            d_cond = 0
            c_condensing = a_cond * (1 / params.dT_approach) + b_cond * (params.T_surface_air_C + 273.15) + c_cond * (params.T_surface_air_C + 273.15) / params.dT_approach + d_cond * (1 / (params.dT_approach + input['dT_range']))
        else:
            raise Exception('CapitalCost_CoolingTower:UnknownCoolingMode', 'Unknown Cooling Mode')

        Q_Ref_BAC = 1e6  # Reference case 1000 kWth (1e6 Wth)
        F_cooling = abs(input['Q_desuperheating']) / (abs(input['Q_desuperheating']) + abs(input['Q_condensing']))
        C_Ref_BAC = Q_Ref_BAC * TDC * (F_cooling * (c_cooling / 1e3) + (1 - F_cooling) * (c_condensing / 1e3))
        C_P = C_Ref_BAC * (abs(input['Q_desuperheating'] + input['Q_condensing']) / Q_Ref_BAC) ** 0.8



        C_P = C_P * (PPI_current / PPI_ref)  # USD

        F_BM = 2.17  # Guthrie Bare-Module Factor for air-fin coolers
        C_BM = C_P * F_BM  # USD

        return C_BM
    
    @staticmethod
    def tank(input, params):
        PPI_ref = PPIs.return_PPI('PPI_Tank', 2013, params)  # 2013 100.00
        PPI_current = PPIs.return_PPI('PPI_Tank', params.costYear, params)  # 2022 232.443

        V = input['V_store']  # m^3
        P = input['P_store'] / 6895  # psi

        P_d = np.exp(0.60608 + 0.91615 * np.log(P) + 0.0015655 * (np.log(P) ** 2))

        S = 15000  # psi
        E = 0.85

        N = 8  # number of vessels
        D = ((V / N) * 4 / 3 / np.pi) ** (1 / 3)  # m
        D = D * 3.28084  # ft
        L = D * 3  # ft

        t_s = (P_d * (D * 12) / (2 * S * E - 1.2 * P_d))  # inches

        rho_V = 0.284  # lb/in^3
        W = np.pi * ((D * 12) + t_s) * ((L * 12) + 0.8 * (D * 12)) * t_s * rho_V  # lb

        C_V = np.exp(5.6336 + 0.4599 * np.log(W) + 0.00582 * (np.log(W)) ** 2)
        C_V = C_V * N

        C_PL = 2275 * D ** 0.2094

        F_BM = 3.05  # Guthrie Bare-Module Factors for vertical pressure vessels
        F_M = 1.0  # carbon steel

        C_P = C_V * F_M + C_PL  # USD
        C_BM = C_P * F_BM  # USD
        C_BM = C_BM * (PPI_current / PPI_ref)  # USD

        return C_BM
    # @staticmethod
    # def cooling_water(input, params):
    #     CEPCI_current = PPIs.return_PPI('CEPCI', params.costYear, params)  # 2022 816

    #     rho = CP.PropsSI('DMASS', 'P', P_pump_in, 'T', T_pump_in)  # kg/m^3
    #     q = m_dot / rho  # m^3/s
        
    #     a = 0.00007 + 2.5e-5 / q
    #     b = 0.003

    #     C_BM_unit = a * CEPCI_current + b  # USD/m^3

    #     C_BM = C_BM_unit * q * 2.88e7  # USD/year
        
    #     return C_BM


    # def refrigeration(input, params):
    #     CEPCI_current = PPIs.return_PPI('CEPCI', params.costYear, params)  # 2022 816

    #     a = 0.5 * (Q_c / 1e3) ** (-0.9) / (T ** 3)
    #     b = 1.1e6 / T ** 5

    #     C_BM_unit = a * CEPCI_current + b  # USD/kJ

    #     C_BM = C_BM_unit * (Q_c / 1e3) * 2.88e7  # USD/year
        
    #     return C_BM


    def TCI(C_TBM):
        C_site = 0.1 * C_TBM
        C_buildings = 0.1 * C_TBM
        C_offsite = 0.05 * C_TBM

        C_TPI = 1.18 * (C_TBM + C_site + C_buildings + C_offsite)

        F_ISF = 1.20
        C_TPI = C_TPI * F_ISF

        C_WC = 1.176 * C_TPI

        C_TCI = C_TPI + C_WC
        
        return C_TCI

