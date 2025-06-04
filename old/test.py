def plot_quivers(x, y, u, v, caption='untitled', verbose=False):
    """plots velocities with PyVista

    x_list - x coordinate array
    y_list - y coordinate array
    us - x velocities at cells
    vs - y velocities at cells  """

    # Create a 2D structured grid object
    grid = pv.StructuredGrid()
  
    grid.dimensions = x.shape[1], x.shape[0], 1     # Set the dimensions of the grid
    
    zip2 = lambda *x: np.array([i for i in zip(*x)])

    # Create the points array with z coordinates as zeros
    z_list = np.zeros_like(x.ravel())
    points = np.stack((x.ravel(), y.ravel(), z_list), axis=-1)
    
    # Set the points of the grid
    grid.points = points
    x_centers = np.array([[0.5*(x[j,i]+x[j,i+1]) for i in range(np.shape(x)[1]-1)] for j in range(np.shape(x)[0]-1)]).ravel()
    y_centers = np.array([[0.5*(y[j,i]+y[j+1,i]) for i in range(np.shape(x)[1]-1)] for j in range(np.shape(x)[0]-1)]).ravel()

    # x_centers = np.array([[0.5(x_list[j,i]+x_list[j,i+1]) for j in range(np.shape(x_list)[0]-1)]])

    speeds = np.pow(np.pow(u, 2) + np.pow(v, 2), 0.5).reshape(np.shape(x)[0]-1, np.shape(y)[1]-1) # plot
    grid.cell_data[caption] = speeds.ravel()
    
    # Visualize the mesh
    plotter = pv.Plotter(off_screen=(not verbose))
    plotter.add_mesh(grid, scalars=caption, cmap='viridis', show_edges=True)
    plotter.add_arrows(zip2(x_centers, y_centers, np.zeros_like(x_centers)), zip2(u, v, np.zeros_like(u)), mag=0.1/np.max(speeds.ravel()))
    plotter.camera_position = 'xy'
    plotter.show_axes()

    # give output
    if not verbose:
        plotter.screenshot(f'{caption}.png')
    elif verbose:
        plotter.show()

    plotter.close()
    del plotter

    return 0