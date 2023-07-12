import numpy as np
import CoolProp.CoolProp as CP
from FrictionFactor import FrictionFactor

def semi_analytic_well(P_f_initial, T_f_initial, T_e_initial, dz_total, dr_total, m_dot, params):
    time_seconds = 3600 * 24 * 365 * params.time_years  # s
    g = 9.81  # m/s^2
    
    # Begin Iterations
    N_dz = 25  # segmentation of the well
    
    # Calculate gradual decreasing inner radius for vertical wells only.
    
    if dr_total == 0 and np.all(np.array(dz_total) < 1):  # injection well; calculating vertical well radii
        wellRadius = np.linspace(params.well_radius, params.side_stream_radius, N_dz)
    elif dr_total == 0 and np.all(np.array(dz_total) > 1):  # production well; calculating vertical well radii
        wellRadius = np.flip(np.linspace(params.well_radius, params.side_stream_radius, N_dz))
    else:  # reservoir; calculating lateral well radius
        wellRadius = np.ones(N_dz) * params.side_stream_radius
    
    dz = dz_total / N_dz  # m
    dr = dr_total / N_dz  # m
    
    dL = np.sqrt(dz**2 + dr**2)  # m
    A_c = np.pi * wellRadius**2  # m^2
    P_c = np.pi * 2 * wellRadius  # m
    
    m = params.n_streams
    
    z = np.zeros((N_dz, m))
    r = np.zeros((N_dz, m))
    L = np.zeros((N_dz, m))
    T_f = np.zeros((N_dz, m))
    T_w = np.zeros((N_dz, m))
    T_e = np.zeros((N_dz, m))
    P = np.zeros((N_dz, m))
    h = np.zeros((N_dz, m))
    rho_fluid = np.zeros((N_dz, m))
    q = np.zeros((N_dz, m))
    Cp_fluid = np.zeros((N_dz, m))
    h_fluid = np.zeros((N_dz, m))
    v = np.zeros((N_dz, m))
    delta_P_loss = np.zeros((N_dz, m))
    
    T_f[0, :] = np.ones(m) * T_f_initial  # C
    T_e[0, :] = np.ones(m) * T_e_initial  # C
    T_w[0, :] = T_f[0, :]  # C
    P[0, :] = np.ones(m) * P_f_initial  # z[0] * 10000; #Pa
    
    h[0, :] = CP.PropsSI('HMASS', 'P', P[0, :], 'T', T_f[0, :] + 273.15, params.fluid)
    rho_fluid[0, :] = CP.PropsSI('DMASS', 'P', P[0, :], 'T', T_f[0, :] + 273.15, params.fluid)
    v[0, :] = m_dot / A_c[0] / rho_fluid[0, :]
    
    # Use Colebrook-white equation for wellbore friction loss.
    # Calculate frictionFactor for first element and assume constant in remainder
    ff = FrictionFactor(wellRadius[0], P[0, :], h[0,:], m_dot, params)

    alpha_rock = params.k_rock / params.rho_rock / params.C_rock  # dim

    # HeatLoss
    # For heat loss calculation
    t_d = alpha_rock * time_seconds / wellRadius ** 2  # dim
    beta = np.zeros(N_dz)
    
    for i in range(N_dz):
        if t_d[i] < 2.8:
            beta[i] = ((np.pi * t_d[i])**-0.5 + 0.5 - 0.25 * (t_d[i] / np.pi)**0.5 + 0.125 * t_d[i])
        else:
            beta[i] = (2 / (np.log(4 * t_d[i]) - 2 * 0.58) - 2 * 0.58 / (np.log(4 * t_d[i]) - 2 * 0.58) ** 2)
    
    for i in range(1, N_dz):
        L[i, :] = L[i - 1, :] + dL
        
        if np.all(np.array(dz_total)) == 0:
            R = dr_total / (2 * np.pi)
            z[i, :] = -(R - R * np.cos(L[i, :] / R)) * np.sin(params.angle * np.pi / 180)
        else:
            z[i, :] = z[i - 1, :] + dz
        
        r[i, :] = r[i - 1, :] + dr
        
        # far-field rock temp
        T_e[i, :] = T_e_initial - z[i, :] * params.dT_dz  # C
        # fluid velocity
        v[i, :] = m_dot / A_c[i - 1] / rho_fluid[i - 1, :]  # m/s
        
        # Calculate Pressure
        delta_P_loss[i, :] = ff * dL / (2 * wellRadius[i - 1]) * rho_fluid[i - 1, :] * v[i, :] ** 2 / 2  # Pa
        rho_fluid[i, :] = CP.PropsSI('DMASS', 'P', P[i - 1, :], 'HMASS', h[i - 1, :], params.fluid)
        P[i, :] = P[i - 1, :] - rho_fluid[i, :] * g * dz - delta_P_loss[i, :]
        
        # Throw exception if below saturation pressure of water at previous temperature
        if params.fluid == 'Water':
            P_sat = CP.PropsSI('P', 'T', T_f[i - 1, :] + 273.15, 'Q', 0, params.fluid)
            if np.any(P[i, :] < P_sat):
                raise Exception('Below saturation pressure of water at {}m !'.format(z[i, :]))
        
        h_noHX = h[i - 1, :] - g * dz
        T_noHX = CP.PropsSI('T', 'P', P[i, :], 'HMASS', h_noHX, params.fluid) - 273.15
        Cp_fluid[i, :] = CP.PropsSI('CPMASS', 'P', P[i, :], 'HMASS', h_noHX, params.fluid)
        

        # Find Fluid Temp
        if not params.useWellboreHeatLoss:
            T_f[i, :] = T_noHX
            h[i, :] = h_noHX
            q[i, :] = 0
        else:
            # See Zhang, Pan, Pruess, Finsterle (2011). A time-convolution
            # approach for modeling heat exchange between a wellbore and
            # surrounding formation. Geothermics 40, 261-266.
            x = dL * P_c[i - 1] * params.k_rock * beta[i - 1] / wellRadius[i - 1]
            y = m_dot * Cp_fluid[i, :]
            
            if np.isinf(np.all(np.array(x))):
                T_f[i, :] = T_e[i, :]
            else:
                T_f[i, :] = (y * T_noHX + x * T_e[i, :]) / (x + y)
            
            q[i, :] = y * (T_noHX - T_f[i, :])
            h[i, :] = CP.PropsSI('HMASS', 'P', P[i, :], 'T', T_f[i, :] + 273.15, params.fluid)
    
    result = {}
    result['State'] = np.concatenate((z[-1, :], P[-1, :] / 1e6, T_f[-1, :], T_e[-1, :]))
    result['stream_depth'] = z[-1, :]
    result['EndTemp'] = T_f[-1, :]
    result['EndPressure'] = P[-1, :]
    result['EndEnthalpy'] = h[-1, :]
    result['EndDensity'] = rho_fluid[-1, :]
    result['Heat'] = -np.sum(q[-1, :])
    result['dP_loss'] = np.sum(delta_P_loss[-1, :])
    result['dL'] = dL
    result['N_dz'] = N_dz
    result['stream_depth'] = z
    result['res_length'] = r
    result['Temp'] = T_f
    result['Pressure'] = P
    result['Enthalpy'] = h
    result['Density'] = rho_fluid
    result['wellRadius'] = wellRadius
    
    return result