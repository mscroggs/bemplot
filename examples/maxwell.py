import bempp.api
import numpy as np
import json

# Solve the (Calderon preconditioned) EFIE
k = 2

def incident_field(x):
    return np.array([0*x[0],0*x[0],np.exp(1j * k * x[0])])

def tangential_trace(x, n, domain_index, result):
    result[:] = np.cross(incident_field(x), n, axis=0)

grid = bempp.api.shapes.sphere(h=0.4)
grid.plot()

e2,e = bempp.api.operators.boundary.maxwell.calderon_electric_field(grid, k)

trace_fun = bempp.api.GridFunction(e.domain, fun=tangential_trace)

from scipy.sparse.linalg import gmres

rhs = e * trace_fun
sol,info = gmres(e2.weak_form(), rhs.coefficients, tol=0.01)
sol = bempp.api.GridFunction(e2.domain, coefficients=sol)

import bemplot

# Export the object itself
pw_c = bempp.api.function_space(grid,"DP",0)
zero = bempp.api.GridFunction(pw_c, coefficients=np.zeros(pw_c.global_dof_count))
bempp.api.export(file_name="max_surface.msh",grid_function=zero)

vm = 4 # the maximum function value. Values larger than this will be clipped.

# define an evaluate. This takes the points as an input and returns the values of the solution at those points
def evaluator(points):
    """Evaluate the off surface solution on the given points."""
    bempp.api.global_parameters.hmat.max_block_size = 2000
    bempp.api.global_parameters.hmat.eps = 1E-5
    slp_pot = bempp.api.operators.potential.maxwell.electric_field(sol.space, points, k)
    scattered_field_data = -slp_pot * sol
    incident_field_data = incident_field(points)
    field_data = scattered_field_data + incident_field_data

    squared_scattered_density = np.real(np.sum(scattered_field_data * scattered_field_data.conj(), axis=0))
    vals = np.real(np.sum(field_data * field_data.conj(), axis=0))
    vals[vals > vm] = vm
    return vals


# Export a slice at x=0 with -6<y<6 and -4<z<8
bemplot.slices("max_slice_x.msh",evaluator,x=[0],ndims=(25,25),extent=[-6,6,-4,8])
# Export a slice at \=0 with -6<x<6 and -6<y<6
bemplot.slices("max_slice_z.msh",evaluator,z=[0],ndims=(25,25),extent=[-6,6,-6,6])

# To make a pretty picture, open max_slice_x.msh, max_slice_z.msh and max_surface.msh with gmsh, then export a screenshot
# For a higher resolution plot, increase ndims above.
