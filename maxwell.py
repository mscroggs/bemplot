from .generic import plot_slice

def efie_evaluator(solution, k, incident=None):
    if incident is None:
        def evaluator(points):
            from bempp.api.operators.potential import maxwell
            elec = maxwell.electric_field(solution.space, points, k)
            plot_me = - elec * solution
            return plot_me
    else:
        def evaluator(points):
            from bempp.api.operators.potential import maxwell
            elec = maxwell.electric_field(solution.space, points, k)
            plot_me = incident(points) - elec * solution
            return plot_me
    return evaluator

def mfie_evaluator(solution, k, incident=None):
    if incident is None:
        def evaluator(points):
            from bempp.api.operators.potential import maxwell
            magn = maxwell.magnetic_field(solution.space, points, k)
            plot_me = - magn * solution
            return plot_me
    else:
        def evaluator(points):
            from bempp.api.operators.potential import maxwell
            magn = maxwell.magnetic_field(solution.space, points, k)
            plot_me = incident(points) - magn * solution
            return plot_me
    return evaluator

def rcfie_evaluator(solution, k, eta, M, incident=None):
    if incident is None:
        def evaluator(points):
            from bempp.api.operators.potential import maxwell
            magn = maxwell.magnetic_field(solution.space, points, k)
            elec = maxwell.electric_field(solution.space, points, k)
            plot_me = - eta * (elec * solution) - (magn * (M * solution))
            return plot_me
    else:
        def evaluator(points):
            from bempp.api.operators.potential import maxwell
            magn = maxwell.magnetic_field(solution.space, points, k)
            elec = maxwell.electric_field(solution.space, points, k)
            plot_me = incident(points) - eta * (elec * solution) - (magn * (M * solution))
            return plot_me
    return evaluator

def cfie_evaluator(solution, k, eta, incident=None):
    if incident is None:
        def evaluator(points):
            from bempp.api.operators.potential import maxwell
            magn = maxwell.magnetic_field(solution.space, points, k)
            elec = maxwell.electric_field(solution.space, points, k)
            plot_me = - eta * (elec * solution) - (magn * solution)
            return plot_me
    else:
        def evaluator(points):
            from bempp.api.operators.potential import maxwell
            magn = maxwell.magnetic_field(solution.space, points, k)
            elec = maxwell.electric_field(solution.space, points, k)
            plot_me = incident(points) - eta * (elec * solution) - (magn * solution)
            return plot_me
    return evaluator
    



def efie_slice(solution, k, plot_type='squared', incident=None, *args, **kwargs):
    plot_slice(efie_evaluator(solution, k, incident), plot_type, *args, **kwargs)

def mfie_slice(solution, k, plot_type='squared', incident=None, *args, **kwargs):
    plot_slice(mfie_evaluator(solution, k, incident), plot_type, *args, **kwargs)

def cfie_slice(solution, k, eta, plot_type='squared', incident=None, *args, **kwargs):
    plot_slice(cfie_evaluator(solution, k, eta, incident), plot_type, *args, **kwargs)

def rcfie_slice(solution, k, eta, M, plot_type='squared', incident=None, *args, **kwargs):
    plot_slice(rcfie_evaluator(solution, k, eta, M, incident), plot_type, *args, **kwargs)

