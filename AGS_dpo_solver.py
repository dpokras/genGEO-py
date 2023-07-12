import time
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import NonlinearConstraint
from SolverFunctions import SolverFunctions as SF
from total_analytic_system_co2 import total_analytic_system_co2

tic = time.time()

def AGS_solver(useMySolver, x0, params):
    if not useMySolver:
        result = total_analytic_system_co2(SF.wrapMyParameters(params, x0))
        x_opt = 0
        print(f"The result is: [{result['W_net']/1e3}] kW")
        print(f"SpCC Power Brownfield: [{result['SpCC_W_brownfield']}] 2019$/W")
        print(f"SpCC Power Greenfield: [{result['SpCC_W_greenfield']}] 2019$/W")
        print(f"SpCC Heat Brownfield: [{result['SpCC_Q_brownfield']}] 2019$/W")
        print(f"SpCC Heat Greenfield: [{result['SpCC_Q_greenfield']}] 2019$/W")
        print(f"Depth: [{params.depth}] m")
        # print(f"Lateral Length: [{x0[1]}] m")
        print(f"Lateral Length: [{params.res_length}] m")
        print(f"no of side streams: [{params.n_streams}]")
        print(f"Mass Flow: [{x0[0]}] kg/s")
        print(f"max speed: [{result['max_speed']}] m/s")
    else: # useMySolver

        myCostFunc = lambda x: SF.costfunc(total_analytic_system_co2(SF.wrapMyParameters(params, x)))
        
        nlc = NonlinearConstraint(lambda x: SF.costfunc(total_analytic_system_co2(SF.wrapMyParameters(params, x))), 0, 3000)

        myOptions = {'maxiter': 20, 'disp': True}

        bounds = [(10, np.inf)]
    
        res = minimize(myCostFunc, x0, method='SLSQP', bounds=bounds, constraints=[nlc], options=myOptions)
        x_opt = res.x

        result = total_analytic_system_co2(SF.wrapMyParameters(params, x_opt))
        print(f"Net Power: [{result['W_net']/1e3}] kW")
        print(f"SpCC Power Brownfield: [{result['SpCC_W_brownfield']}] 2019$/W")
        print(f"SpCC Heat Brownfield: [{result['SpCC_Q_brownfield']}] 2019$/W")
        print(f"SpCC Power Greenfield: [{result['SpCC_W_greenfield']}] 2019$/W")
        print(f"SpCC Heat Greenfield: [{result['SpCC_Q_greenfield']}] 2019$/W")
        print(f"Depth: [{params.depth}] m")
        print(f"Lateral Length: [{params.res_length}] m")
        # print(f"Lateral Length: [{x_opt[1]}] m")
        print(f"no of side streams: [{params.n_streams}]")
        print(f"Mass Flow: [{x_opt[0]}] kg/s")
        print(f"max speed: [{result['max_speed']} m/s")
        print(f"entopy into turbine: [{result['s_turb_in']}] J/kgK")
    toc = time.time()

    print(toc - tic)
    return result, x_opt