def slices(file_name, evaluator, extent=(-5,5,-5,5), x=None,y=None,z=None, ndims=(300,300)):
    if x is not None:
        mode = "yz"
        ny,nz = ndims
        ps = x
        axis="x"
    elif y is not None:
        mode = "xz"
        nx,nz = ndims
        ps = y
        axis="y"
    elif z is not None:
        mode = "xy"
        nx,ny = ndims
        ps = z
        axis="z"
    else:
        raise ValueError("One of x, y and z must be set")
    import bempp.api
    import os

    fname, extension = os.path.splitext(file_name)

    ll = (extent[0],extent[2])
    ur = (extent[1],extent[3])

    node_offset = 1
    element_offset = 1

    for p in ps:
        grid = bempp.api.structured_grid(ll, ur, ndims, axis=mode, offset=p)
        nnodes = grid.leaf_view.entity_count(2)
        nelements = grid.leaf_view.entity_count(0)
        space = bempp.api.function_space(grid, "P", 1, domains=[0], closed=True)
        points = space.global_dof_interpolation_points
        vals = evaluator(points)
        output_fun = bempp.api.GridFunction(space, coefficients=vals)
        bempp.api.export(file_name=fname + "_" + axis + "=" + str(p) + extension, grid_function=output_fun,
                         data_type='node',
                         vertex_index_to_file_key_map=range(
                             node_offset, node_offset + nnodes),
                         element_index_to_file_key_map=range(
                             element_offset, element_offset + nelements))
        node_offset += nnodes
        element_offset += nelements
