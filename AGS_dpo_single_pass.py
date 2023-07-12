from PPIs import PPIs
from SimulationParameters import SimulationParameters
import pandas as pd
from AGS_dpo_solver import AGS_solver as solver

# Define System Parameters
params = SimulationParameters()

params.time_years = 30
params.wellCostType = 'Baseline'
params.config = 1

params.well_radius = 0.4667 / 2  # m ~18 3/8" inner diameter
params.side_stream_radius = 0.2032 / 2  # m ~8" inner diameter    
params.S_ratio = 0.95  # S_ratio = m_S1/m_total

if params.config == 1 or params.config == 2:
    params.S_ratio = 1  # S_ratio = m_S1/m_total
    
params.find_opt_S_ratio_toggle = 1
params.coolingMode = 'Wet'

# Solver toggle
useMySolver = 0
x0 = [90]

# Calculate
params.n_streams = 11
params.res_length = 7000

result, x_opt = solver(useMySolver, x0, params)

if result['s_turb_in'] <= (1.65e3 * 0.99):
    print('constraint NOT satisfied')
    con = 0
else:
    print('constraint satisfied')
    con = 1

A = [result['W_net'] / 1e3,
     result['SpCC_W_brownfield'],
     result['SpCC_Q_brownfield'],
     result['SpCC_W_greenfield'],
     result['SpCC_Q_greenfield'],
     params.depth,
     params.res_length,
     params.n_streams,
     x0,
     result['max_speed'],
     result['CapitalCost']['C_brownfield'],
     result['CapitalCost']['C_greenfield'],
     con]

T = pd.DataFrame([A], columns=['W_net',
                               'SpCC_W_Br',
                               'SpCC_Q_Br',
                               'SpCC_W_Gr',
                               'SpCC_Q_Gr',
                               'depth',
                               'res_length',
                               'n_streams',
                               'm_dot',
                               'max_speed',
                               'C_Br',
                               'C_Gr',
                               'con satisfied'])