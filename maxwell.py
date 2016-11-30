def efie_slice(solution, k, plot_type='squared', incident=None, x=None, y=None, z=None, n=(151,151), extent=(-5,5,-5,5), cmap='coolwarm', filename=None, title=None, show=None, **kwargs):
    from bempp.api.operators.potential import maxwell
    import numpy as np

    if show is None:
        show = filename==None

    if x is not None and y is None and z is None:
        x_p, y_p, z_p = np.mgrid[x:x:1j, extent[0]:extent[1]:n[0]*1j, extent[2]:extent[3]:n[1] * 1j]
    elif x is None and y is not None and z is None:
        x_p, y_p, z_p = np.mgrid[extent[0]:extent[1]:n[0]*1j, y:y:1j, extent[2]:extent[3]:n[1] * 1j]
    elif x is None and y is None and z is not None:
        x_p, y_p, z_p = np.mgrid[extent[0]:extent[1]:n[0]*1j, extent[2]:extent[3]:n[1] * 1j, z:z:1j]
    else:
        raise TypeError("Exactly one of x, y and z must be set.")
    points = np.vstack((x_p.ravel(), y_p.ravel(), z_p.ravel()))

    elec = maxwell.electric_field(solution.space, points, k)

    if incident is None:
        plot_me = - elec * solution
    else:
        plot_me = incident(points) - elec * solution

    import matplotlib
    from matplotlib import pyplot as plt

    if plot_type == 'squared':
        plot_me = np.real(np.sum(plot_me * plot_me.conj(), axis=0))
    elif plot_type == 'x':
        plot_me = np.real(plot_me[0,:])
    elif plot_type == 'y':
        plot_me = np.real(plot_me[1,:])
    elif plot_type == 'z':
        plot_me = np.real(plot_me[2,:])
    else:
        raise ValueError("plot_type value invalid")

    plt.imshow(plot_me.reshape(n).T,
               cmap=cmap, origin='lower',
               extent=extent, **kwargs)
    plt.colorbar()
    if title is None:
        if x is not None:
            title = "Plot at x="+str(x)
        if y is not None:
            title = "Plot at y="+str(y)
        if z is not None:
            title = "Plot at z="+str(z)
    plt.title(title)
    if x is None:
        plt.xlabel("x")
        if y is None:
            plt.ylabel("y")
        else:
            plt.ylabel("z")
    else:
        plt.xlabel("y")
        plt.ylabel("z")
    if filename is not None:
        plt.savefig(filename)
    if show:
        plt.show()

    plt.clf()

