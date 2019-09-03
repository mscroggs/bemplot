import bempp.api
import numpy as np

# Solve a Helmholtz problem
k = 6.

def incident(x):
    return np.exp(1j * k * (x[0]+x[1]+x[2])/np.sqrt(3))

def trace(x, n, domain_index, result):
    result[:] = incident(x)

grid = bempp.api.shapes.cube(h=0.4)
space = bempp.api.function_space(grid, "P", 1)

id = bempp.api.operators.boundary.sparse.identity(space, space, space)
V = bempp.api.operators.boundary.helmholtz.single_layer(space, space, space, k)
K = bempp.api.operators.boundary.helmholtz.double_layer(space, space, space, k)

trace_fun = bempp.api.GridFunction(space, fun=trace)

sol,info = bempp.api.linalg.gmres(V, -(-0.5*id+K)*trace_fun)


import bemplot

# Export the object itself
pw_c = bempp.api.function_space(grid,"DP",0)
zero = bempp.api.GridFunction(pw_c, coefficients=np.zeros(pw_c.global_dof_count))
bempp.api.export(file_name="helm_surface.msh",grid_function=zero)

vm = 4 # the maximum function value. Values larger than this will be clipped.

# define an evaluate. This takes the points as an input and returns the values of the solution at those points
def evaluator(points):
    """Evaluate the off surface solution on the given points."""
    bempp.api.global_parameters.hmat.max_block_size = 2000
    bempp.api.global_parameters.hmat.eps = 1E-5
    slp_pot = bempp.api.operators.potential.helmholtz.single_layer(sol.space, points, k)
    scattered_field_data = - slp_pot * sol
    incident_field_data = incident(points)
    field_data = scattered_field_data + incident_field_data

    vals = np.real(np.sum(field_data * field_data.conj(), axis=0))
    vals[vals > vm] = vm
    return vals


# Export a slices at x=0 and x=0.5 with -6<y<6 and -4<z<8
bemplot.slices("helm_slice_x.msh",evaluator,x=[0,0.5],ndims=(25,25),extent=[-6,6,-4,8])
# Export a slice at z=0 with -6<x<6 and -6<y<6
bemplot.slices("helm_slice_z.msh",evaluator,z=[0.5],ndims=(25,25),extent=[-6,6,-6,6])

# To make a pretty picture, open helm_slice_x=0.msh, helm_slice_z=0.5.msh and helm_surface.msh with gmsh, then export a screenshot
# For a higher resolution plot, increase ndims above.
