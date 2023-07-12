import numpy as np

def LCOE_Simple(CapitalCost, Capacity_W, Capacity_Q, params):
    # Capacity input in Watts
    # Capital cost input in $
    # LCOE output in $/W-h
    # SpCC output in $/W
    
    # Check heats & powers
    if Capacity_W <= 0 or CapitalCost <= 0:
        return np.nan, np.nan
    
    if params.F_OM < 0:
        raise Exception('Negative O&M Fraction')
    
    if params.discountRate < 0:
        raise Exception('Negative Discount Rate')
    
    if params.Lifetime <= 0:
        raise Exception('Negative Lifetime')
    
    if params.CapacityFactor <= 0 or params.CapacityFactor > 1:
        raise Exception('Bad Capacity Factor')
    
    CRF = params.discountRate * (1 + params.discountRate)**params.Lifetime / ((1 + params.discountRate)**params.Lifetime - 1)
    
    LCOE = CapitalCost * (CRF + params.F_OM) / (Capacity_W * params.CapacityFactor * 8760)
    if Capacity_W != 0:
        SpCC_W = CapitalCost / Capacity_W
    else:
        SpCC_W = 0
    if Capacity_Q != 0:
        SpCC_Q = CapitalCost / Capacity_Q
    else:
        SpCC_Q = 0
        
    result = {}
    result['LCOE'] = LCOE
    result['SpCC_W'] = SpCC_W
    result['SpCC_Q'] = SpCC_Q
    
    return result