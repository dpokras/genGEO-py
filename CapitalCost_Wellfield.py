from CapitalCost_Well import CapitalCost_Well
from PPIs import PPIs

def CapitalCost_Wellfield(params):
    PPI_Permit = PPIs.return_PPI('PPI_Permit', params.costYear, params)  # 2022
    PPI_OG_s = PPIs.return_PPI('PPI_O&G-s', params.costYear, params)  # 2022
    X_IC_wf = 1.05
    X_PC_wf = 1.15

    WellFieldLength_km = 1  # Assuming it's a square wellfield
    G_injection = WellFieldLength_km ** 2
    G_monitoringWells = G_injection

    PPI_Permit_ref = PPIs.return_PPI('PPI_Permit', 2003, params)  # 2003
    C_permitting = X_IC_wf * X_PC_wf * PPI_Permit / PPI_Permit_ref * 665700

    L_wellZone = WellFieldLength_km * 1000
    A_CO2_AOR = L_wellZone ** 2
    A_CO2_AMA = (L_wellZone + 1600) ** 2
    dC_permitting_CO2 = X_IC_wf * X_PC_wf * PPI_Permit / PPI_Permit_ref * 45000 * (1 / 1e6) * A_CO2_AMA

    PPI_OG_s_ref = PPIs.return_PPI('PPI_O&G-s', 2003, params)  # 2003
    C_monitoringWell = CapitalCost_Well(params.depth, 0, params.monitoringWellDiameter / 2, params)
    C_monitoringWells_CO2 = G_monitoringWells * C_monitoringWell
    C_surfaceMonitoring_CO2 = X_IC_wf * X_PC_wf * PPI_OG_s / PPI_OG_s_ref * 138000 * (1 / 1e6) * A_CO2_AMA
    dC_monitoring_CO2 = C_monitoringWells_CO2 + C_surfaceMonitoring_CO2

    C_wellfield = C_permitting

    return C_wellfield
