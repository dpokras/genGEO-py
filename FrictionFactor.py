import numpy as np
import CoolProp.CoolProp as CP

def FrictionFactor(wellRadius, P, h, m_dot, params):
    rho_fluid = CP.PropsSI('DMASS', 'P', P, 'HMASS', h, params.fluid)
    mu = CP.PropsSI('VISCOSITY', 'P', P, 'HMASS', h, params.fluid)

    A_c = np.pi * wellRadius**2
    V = m_dot / A_c / rho_fluid

    # Relative Roughness
    rr = params.epsilon / (2 * wellRadius)

    # Reynolds Number
    Re = rho_fluid * V * (2 * wellRadius) / mu

    # Use Haaland (1983) Approximation
    c1 = -1.8 * np.log10((6.9 / Re) + (rr / 3.7)**1.11)
    ff = (1 / c1)**2

    return ff

