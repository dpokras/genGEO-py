from PPIs import PPIs

def CapitalCost_Well(welldepth, wellLength, wellradius, params):
    PPI_current = PPIs.return_PPI('PPI_O&G', params.costYear, params)  # 2022 2.195
    
    X_IC_well = 1.05
    X_PC_well = 1.15
    
    if welldepth > 0:
        if params.wellCostType == 'Ideal':
            # revised GETEM cost correlations 2022 large diameter -
            # 12.25in vertical well)
            PPI_ref = PPIs.return_PPI('PPI_O&G', 2003, params)  # 2003 153.492
            C_well = X_IC_well * X_PC_well * PPI_current / PPI_ref * (-62.2 * welldepth + 1290 * (2 * wellradius) * welldepth + 275300)
            # disp('Caution: Using Ideal Well Relations')
        elif params.wellCostType == 'Baseline':
            PPI_ref = PPIs.return_PPI('PPI_O&G', 2003, params)  # 2003 153.492
            C_well = X_IC_well * X_PC_well * PPI_current / PPI_ref * (0.105 * welldepth ** 2 + 1776 * (2 * wellradius) * welldepth + 275300)
        else:
            raise ValueError('Unknown Well Type')
        dC_well_CO2 = X_IC_well * X_PC_well * PPI_current / PPI_ref * (265 * (2 * wellradius) * welldepth + 133 * welldepth)
    else:
        # horizontal well (revised GETEM correlation 2022 small diameter -
        # 8.5 in horizontal well)
        PPI_ref = PPIs.return_PPI('PPI_O&G', 2003, params)  # 2003 153.492
        # horizontal well cost modifier = 15%
        X_HWell = 1.15
        C_well = X_IC_well * X_PC_well * PPI_current / PPI_ref * X_HWell * (0.105 * wellLength ** 2 + 1776 * (2 * wellradius) * wellLength + 275300)
        dC_well_CO2 = X_IC_well * X_PC_well * PPI_current / PPI_ref * (265 * (2 * wellradius) * wellLength + 133 * wellLength)
    
    if params.fluid == 'Water':
        WellCost = C_well / params.SuccessRate_well
    elif params.fluid == 'CO2':
        WellCost = (C_well + dC_well_CO2) / params.SuccessRate_well
    else:
        raise NotImplementedError('Not Implemented')
    
    return WellCost
