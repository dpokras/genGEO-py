class SolverFunctions:
    def wrapMyParameters(params, x):
        params.m_dot = x[0]
        # params.res_length = x[1]
        return params

    def costfunc(result):
        cost = result['SpCC_W_greenfield']
        return cost

    def contrfunc(result):
        c = 1.65e3 - result['s_turb_in']
        # c = 100 - result['T_prod_surface']
        ceq = []
        # ceq = result['W_net'] / 1e3 - 1.2e3
        return c, ceq