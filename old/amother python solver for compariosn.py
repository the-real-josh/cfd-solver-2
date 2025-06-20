import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
import pyvista as pv
from numpy import float64
import os 

epsilon = 1e-7
savedir = r'results\\'

def find_arg_nearest(array,value):
    """inputs: the array you want to search, value that you want to seach for
        outputs: the index of the value closest to the value you searched for"""
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or np.fabs(value - array[idx-1]) < np.fabs(value - array[idx])):
        return idx-1
    else:
        return idx


def mesh_plot(x_list, y_list, title, return_squishedness=False):
    """Pass data in as x_list and y_list as a list of lists in the shape of the mesh"""
        # fancy plot with skewness
    squishedness = 0
    n = 0
    i_ubound = len(x_list)
    j_ubound = len(x_list[0])
    for i in range(i_ubound):
        for j in range(j_ubound):
            try:
                if (i==i_ubound and j==j_ubound) or (i==i_ubound-1 and j==j_ubound-1):
                    pass
                elif j==j_ubound-1:
                    plt.plot((x_list[i, j], x_list[i+1, j]), (y_list[i, j], y_list[i+1, j]), color='black')
                elif i==i_ubound-1:
                    plt.plot((x_list[i, j], x_list[i, j+1]), (y_list[i, j], y_list[i, j+1]), color='black')
                else:
                    plt.plot((x_list[i, j], x_list[i+1, j]), (y_list[i, j], y_list[i+1, j]), color='black')
                    plt.plot((x_list[i, j], x_list[i, j+1]), (y_list[i, j], y_list[i, j+1]), color='black')
                    horiz = np.array([x_list[i+1, j], y_list[i+1, j]]) - np.array([x_list[i, j], y_list[i, j]])
                    vert = np.array([x_list[i, j+1], y_list[i, j+1]]) - np.array([x_list[i, j], y_list[i, j]])
                    if not 0:
                        squishedness += np.abs(np.dot(horiz, vert) / (np.linalg.norm(horiz)*np.linalg.norm(vert)))
                        n += 1
            except IndexError:
                print(f'grapher - index error @ji={j,i} for shape of {np.shape(x_list)}')
                pass
    plt.title(title)
    plt.scatter((2, 2, 3, 3), (0, 1, 0, 1))
    plt.axis('equal')
    plt.show()
    if return_squishedness:
        return squishedness/(i_ubound*j_ubound)


def mesh(mesh_shape, min_residual=1e-6):
    """# Resolution:
        - resolution of under 4 causes everything to completely break
        - resolution is nodes per unit
        - change the resolution to prove mesh-independence of solution

        # min residual:
         - the average amount that is changed by the iterator
         - default is 1e-6
         """
    print(f'\n{mesh_shape} mesh not found - generating a {mesh_shape}...')

    # DEFINE FUNCTIONS FOR THE TOP AND BOTTOM OF CV
    def y_lower(x):
        # return 0
        if 0 <= x <= 2:
            return 0.0
        elif 2 < x <= 3:
            return 0.12*np.sin((x-2)*3.14159)
        elif 3 < x <= 5:
            return 0.0
        else:
            return np.nan
        
    def y_upper(x):
        # return 1
        if 0 <= x <= 2:
            return 1.0
        elif 2 < x <= 3:
            return 1.0-0.12*np.sin((x-2)*3.14159)
        elif 3 < x <= 5:
            return 1.0
        else:
            return np.nan

    # DEFINE BOUNDS
    x_l = 0; x_r = 5 # left and right bound of x
    J_max = mesh_shape[0] # number of nodes in the y direction
    I_max = mesh_shape[1] # for debuggin'


    # generate xi-eta grid

    # xi (x)
    xi_list = np.zeros((J_max, I_max), dtype=float) # create J_max rows of I_max entries each
    for ari in range(J_max): # for every row in the J_max steps
        for i in range(I_max): # generate a row with I_max steps
            xi_list[ari, i] = (i) / (I_max-1)

    # eta (y)
    eta_list = np.zeros((J_max, I_max), dtype=float) # create __ rows of ___ entries each
    for i in range(I_max): # for every row,
        for ari in range(J_max): # generae a row.
            eta_list[ari, i] = (ari) / (J_max-1)
    
    # test plot
    # plt.title(f'Algebraic Grid in ξ-η Space'); plt.plot(xi_list, eta_list, marker='.', color='k', linestyle='none'); plt.show()

    x_list = np.zeros((J_max, I_max), dtype=float); y_list = x_list.copy() # pre-allocate memory. Remember - zeroes works by generating __ rows of ___ length
    # print(f'shape of x list: {np.shape(x_list)}')
    for i in range(I_max):
        for ari in range(J_max):
            try:
                x_list[ari, i] = x_l + xi_list[ari, i] * (x_r - x_l) # linearly scale x up
            except:
                print(f'bad x, {[i, ari]} out of bounds for either x list: {np.shape(x_list)} or xi_list: {np.shape(xi_list)}')
                print(x_list)
                exit()
            try:
                x = x_list[0, i] # COULD take advantage of the fact that the domain on the x scales from 0 to 1 and just say xi, but eh
                y_list[ari, i] = y_lower(x) + eta_list[ari, i] * (y_upper(x) - y_lower(x))
            except:
                print(f'bad y')
                print(y_list)
                exit()

    # initial algebraic grid plot
    squishedness = mesh_plot(x_list, y_list, 'Algebraic grid', return_squishedness=True)
    print(f'Algebraic Average Skewness - (lower is better): {squishedness}')
    plt.show()

    # preliminary calcs
    delta_xi = 1/(I_max);   delta_eta = 1/(J_max) # this is because xi and eta are defined to be on the interval {0<x<1}. In python if your list has max index of 3, then you have 2 divisions

    old_xlist = x_list.copy()
    old_ylist = y_list.copy()


    # break out of the loop either when residuals hit the floor or iterations hit the ceiling
    max_iterations = 100*mesh_shape[0] # iteration ceiling
    # min_residual = 1e-6 # residual floor

    aggregate_residuals = np.zeros((max_iterations+1))
    print(f'Beginning a maximum of {max_iterations} iterations.\nCurrent iteration:  ', end='')

    iterations = 0
    while iterations < max_iterations:
        iterations += 1

        # fancy schmancy status
        print('\b'*len(str(iterations-1)), end='', flush=True)
        print(iterations, end='', flush=True)


        corner_pts = [find_arg_nearest(x_list[0], 2.0), find_arg_nearest(x_list[0], 3.0)] # define the corner points and the bool function (for simplicity)
        def corner_point(i,j):
            if j in [0, J_max-1] and i in corner_pts:
                return True
            else:
                return False

        # sets the iterations for i,j
        iter_set = np.array([[(i, j) for i in range(I_max)] for j in range(J_max)]).reshape(I_max*J_max, 2)

        iter_residual = 0
        for i,j in iter_set:
            if (corner_point(i,j)):
                if j == 0 and i < I_max/2:
                    new_x_val = 2.0
                    new_y_val = 0.0
                    y_list[j,i] = new_y_val
                elif j == J_max - 1 and i < I_max/2:
                    new_x_val = 2.0
                    new_y_val = 1.0
                    y_list[j,i] = new_y_val
                elif j == 0 and i > I_max/2:
                    new_x_val = 3.0
                    new_y_val = 0.0
                elif j == J_max-1 and i > I_max/2:
                    new_x_val = 3.0
                    new_y_val = 1.0
            elif j in [0, J_max - 1]: # for top and bottom
                if x_list[j,i] < 2 or x_list[j,i] > 3:
                    try: 
                        new_x_val = x_list[j+1,i]
                    except:
                        new_x_val = x_list[j-1,i]
                    new_x_val = x_list[j,i] # override
                if x_list[j,i] < 3 and x_list[j,i] > 2:
                    try:
                        new_x_val = x_list[j+2, i] + ((y_lower(x_list[j,i-1]) - y_lower(x_list[j,i+1]))/(x_list[j,i-1] - x_list[j,i+1])) * delta_xi * x_r 
                        new_y_val = y_lower(new_x_val)
                        iter_residual += abs(y_list[j,i] - new_y_val)
                        y_list[j,i] = new_y_val
                    except:
                        new_x_val = x_list[j-2, i] + ((y_lower(x_list[j,i-1]) - y_lower(x_list[j,i+1]))/(x_list[j,i-1] - x_list[j,i+1])) * delta_xi * x_r 
                        new_y_val = y_upper(new_x_val)
                        iter_residual += abs(y_list[j,i] - new_y_val)
                        y_list[j,i] = new_y_val
                iter_residual += abs(x_list[j,i] - new_x_val)
                x_list[j,i] = new_x_val
            elif i in [0, I_max - 1]:
                if i == 0:
                    new_y_val = y_list[j,i + 1]
                    iter_residual += abs(y_list[j,i] - new_y_val)
                elif i == I_max - 1:
                    new_y_val = y_list[j,i - 1]
                iter_residual += abs(y_list[j,i] - new_y_val)
                y_list[j,i] = new_y_val
            
            else: # general inside equations
                alpha = 1/(4*delta_eta**2)  * ((x_list[j+1,i] - x_list[j-1,i])**2 + (y_list[j+1,i] - y_list[j-1,i])**2) # see cizmas p.33
                beta = -1/(4*delta_eta*delta_xi) * ((x_list[j,i+1] - x_list[j,i-1])*(x_list[j+1,i] - x_list[j-1,i])   +   (y_list[j,i+1] - y_list[j,i-1])*(y_list[j+1,i] - y_list[j-1,i]))
                gamma = 1/(4*delta_xi**2) * ((x_list[j,i+1] - x_list[j,i-1])**2 + (y_list[j,i+1] - y_list[j,i-1])**2)

                k1 = alpha / delta_xi**2
                k2 = -beta / (2*delta_xi*delta_eta)
                k3 = gamma / delta_eta**2


                new_x_val = (0.5/(k1+k3)) * (k2*x_list[j-1,i-1] + k3*x_list[j-1,i] -k2*x_list[j-1,i+1]+ k1*x_list[j,i-1]+ k1*x_list[j,i+1] -k2*x_list[j+1,i-1] + k3*x_list[j+1,i] + k2*x_list[j+1,i+1])
                x_list[j,i] = new_x_val
                iter_residual += abs(x_list[j,i] - new_x_val)

                new_y_val = (0.5/(k1+k3)) * (k2*y_list[j-1,i-1] + k3*y_list[j-1,i] -k2*y_list[j-1,i+1]+ k1*y_list[j,i-1]+ k1*y_list[j,i+1] -k2*y_list[j+1,i-1] + k3*y_list[j+1,i] + k2*y_list[j+1,i+1])
                iter_residual += abs(y_list[j,i] - new_y_val)
                y_list[j,i] = new_y_val
                

        aggregate_residuals[iterations] = iter_residual/(I_max * J_max)
        if aggregate_residuals[iterations] < min_residual:
            print(f'\nResidual<{min_residual} - convergence detected!')
            break
    else:
        print(f'\niteration ceiling reached')

    # plot them
    plt.title(f'Test plot of Elliptic Grid Vertices in XY at 0 and {iterations} iterations'); plt.plot(x_list, y_list, marker='.', color='k', linestyle='none'); plt.axis('equal')
    plt.plot(old_xlist, old_ylist, marker='.', color='r', linestyle='none'); plt.axis('equal'); 
    plt.show()

    squishedness = mesh_plot(x_list, y_list, 'Elliptic grid', return_squishedness=True)
    print(f'Elliptic Grid Skewness - lower is better: {squishedness}')

    # plot residuals
    aggregate_residuals = aggregate_residuals[1:iterations]
    plt.plot(np.linspace(0, iterations, num=len(aggregate_residuals)), aggregate_residuals)
    plt.title(f'Residuals')
    plt.ylabel(f'Residual mean value on iteration')
    plt.xlabel(f'Iteration number')
    plt.yscale('log')
    plt.show()
    plt.yscale('linear')

    print(f'\nmesh generation complete.')
    return pd.DataFrame({f'{mesh_shape}-x':x_list.ravel(),
                            f'{mesh_shape}-y':y_list.ravel()})


def trim_dataframe(df):
    if np.isnan(df.iloc[-1]):
        # Find the last non-blank entry index for each column
        last_good = np.argmax(df.apply(lambda col: np.isnan(col)).to_numpy())
        trimmed_df = df[:last_good]
        return trimmed_df
    else:
        return df

def plot_quivers(x_list, y_list, us, vs, caption='untitled', verbose=False):
    """plots velocities with PyVista

    x_list - x coordinate array
    y_list - y coordinate array
    us - x velocities at cells
    vs - y velocities at cells  """

    # Create a 2D structured grid object
    grid = pv.StructuredGrid()
  
    grid.dimensions = x_list.shape[1], x_list.shape[0], 1     # Set the dimensions of the grid
    
    zip2 = lambda *x: np.array([i for i in zip(*x)])

    # Create the points array with z coordinates as zeros
    z_list = np.zeros_like(x_list.ravel())
    points = np.stack((x_list.ravel(), y_list.ravel(), z_list), axis=-1)
    
    # Set the points of the grid
    grid.points = points
    x_centers = np.array([[0.5*(x_list[j,i]+x_list[j,i+1]) for i in range(np.shape(x_list)[1]-1)] for j in range(np.shape(x_list)[0]-1)]).ravel()
    y_centers = np.array([[0.5*(y_list[j,i]+y_list[j+1,i]) for i in range(np.shape(x_list)[1]-1)] for j in range(np.shape(x_list)[0]-1)]).ravel()

    # x_centers = np.array([[0.5(x_list[j,i]+x_list[j,i+1]) for j in range(np.shape(x_list)[0]-1)]])

    speeds = np.pow(np.pow(us, 2) + np.pow(vs, 2), 0.5).reshape(np.shape(x_list)[0]-1, np.shape(y_list)[1]-1) # plot
    grid.cell_data[caption] = speeds.ravel()
    
    # Visualize the mesh
    plotter = pv.Plotter(off_screen=(not verbose))
    plotter.add_mesh(grid, scalars=caption, cmap='viridis', show_edges=True)
    plotter.add_arrows(zip2(x_centers, y_centers, np.zeros_like(x_centers)), zip2(us, vs, np.zeros_like(us)), mag=0.1/np.max(speeds.ravel()))
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

def plot_scalars(x_list, y_list, scalar_data, caption='untitled', verbose=False):
    """plots with PyVista"""
    # use pyvista (https://docs.pyvista.org/user-guide/)

    # Create a 2D structured grid object
    grid = pv.StructuredGrid()
 
    grid.dimensions = x_list.shape[1], x_list.shape[0], 1     # Set the dimensions of the grid
    
    pv.set_plot_theme("document")
    grid.dimensions = x_list.shape[1], x_list.shape[0], 1     # Set the dimensions of the grid
    
    # Create the points array with z coordinates as zeros
    z_list = np.zeros_like(x_list.flatten())
    points = np.stack((x_list.flatten(), y_list.flatten(), z_list), axis=-1)
    
    # Set the points of the grid
    grid.points = points
    grid.cell_data[caption] = scalar_data.ravel()

    # only put edges when the cells are reasonably large
    if len(x_list) > 30:
        edges = False
    else:
        edges = True

    plotter = pv.Plotter(off_screen=(not verbose))
    plotter.add_mesh(grid, scalars=caption, cmap='viridis', show_edges=edges)
    plotter.camera_position = 'xy'
    plotter.show_axes()

    # give output
    if not verbose:
        plotter.screenshot(f'{savedir}{caption}.png')
    elif verbose:
        plotter.show()

    plotter.close()
    del plotter

    return 0


def run(x_list, y_list, mach):

    # append "ghost cells"
    print(f'initializing ghost cells: cell size: {np.shape(x_list)}-->', end=' ')
    delta_x = np.average(np.diff(x_list[0, :]))
    delta_y = np.average(np.diff(y_list[:, 0].ravel()))

    # append in the top and bottom
    x_list = np.concatenate(([x_list[0, :]], [x_list[0, :]], x_list[:, :], [x_list[-1, :]], [x_list[-1, :]]), axis=0)
    y_list = np.concatenate(([y_list[0, :] - 2*delta_y], [y_list[0, :] - delta_y], y_list[:, :], [y_list[-1, :] + delta_y], [y_list[-1, :] + 2*delta_y]), axis=0)

    # append in the left and the right
    x_list = np.concatenate((x_list[:, 0:1] - 2*delta_x, x_list[:, 0:1] - delta_x, x_list[:], x_list[:,  -1:] + delta_x, x_list[:,  -1:] + 2*delta_x), axis=1)
    y_list = np.concatenate((y_list[:, 0:1], y_list[:, 0:1], y_list[:], y_list[:,  -1:], y_list[:,  -1:]), axis=1)

    print(f'{np.shape(x_list)}')
    
    # initial gas properties - only need 2 properties to fix a simple state
    rho_infty = 1.212 # (kg/m3)
    T_infty = 293 # if you took very warm air then expanded it to current mach

    # fixed properties 
    gamma = 1.4 # ratio of specific heats
    R = 287.052874 # gas constant for air, J/(kg K)
    cp = gamma*R / (gamma-1)
    cv = cp - R
    a = np.sqrt(gamma * R * T_infty) # speed of sound (m/s)
    # E_inlet = 1.0/(rho * (gamma-1)) + 0.5 * mach**2 
    # E_inlet = 1.0/(gamma * (gamma-1)) + 0.5 * mach**2 # cizmas in the textbook 
    V_inlet = a*mach
    E_inlet = T_infty*cv + 1/2 * (V_inlet)**2 # total energy at the inlet, J/kg
    p_infty = (T_infty*cv)*rho_infty * (gamma - 1) # use static energy (not stagnation)
    print(f'initial speed of sound: {a} m/s\n \
            initial temperature: {T_infty} K\n\
            initial density: {rho_infty} kg/3\n\
            initial energy: {E_inlet} J/kg')

    # 2ND ORDER DISSIPATION COEFFICIENT
    # cizmas and TA say 0
    # cizmas notes say 0-0.5
    # Jameson says 0.25
    nu_2 = 0.25

    # 4TH ORDER DISSIPATION COEFFICIENT
    #nu4 - cizmas says 0.001 or 0.01, Jameson says 1/256
    # nu_4 = 0.00390625 # or 0.01 or 0.001
    # TAs say 0.01 
    # 0.02 works for Mach 0.3 at 11x51 grid size
    nu_4 = 0.00390625

    # any higher and it risks blowing up. Should be able to put CFL = 1 with no issue.
    CFL = 0.2
    alphas = [1/4, 1/3, 1/2, 1.0]

    # ======= INITIALIZE THE STATE ============

    # BEGIN WITH STATE VECTOR Q
    try:
        # try to recover previous work
        filename = f'{savedir}j_max={len(x_list)-1},i_max={len(x_list[0])-1},M={round(mach, 1)}.csv'
        df = pd.read_csv(filename)
        cell_data = df['data'].dropna().to_numpy().reshape(len(x_list)-1, len(x_list[0])-1, 4)
        ag_res_list = df['convergence'].dropna().to_list()
        it_number = len(ag_res_list)
        print(f'Found file {filename}\nContinuing from iteration {it_number}')
    except FileNotFoundError:
        print(f'Could not find file with correct specifications in cwd. Starting from scratch.\n')
        it_number = 0 # n is the iteration number
        def theta_cell(j,i):
            try:
                return np.tan((0.5*(y_list[j,i+1] - y_list[j,i-1] + y_list[j+1,i+1] - y_list[j+1,i-1]))/(x_list[j,i+1] - x_list[j,i-1]))
            except IndexError:
                return 0.0
        # uniform following cell's geometry
        cell_data = np.array([[np.array([rho_infty, rho_infty*V_inlet*np.cos(theta_cell(j,i)), rho_infty*V_inlet*np.sin(theta_cell(j,i)), rho_infty*E_inlet], dtype=float) for i in range(len(x_list[0])-1)] for j in range(len(x_list)-1)]) # initialize cell data - start with uniform rightward flow
        
        # uniform rightward
        # cell_data = np.array([[np.array([rho_infty, rho_infty*V_inlet, 0.0, rho_infty*E_inlet], dtype=float) for i in range(len(x_list[0])-1)] for j in range(len(x_list)-1)])
        ag_res_list = [] # list of all the residuals ever
        
    tot_mass_list = []
    avg_u_list = []
    avg_v_list = []
    tot_energy_list = []
        

    # pre-allocate tools
    new_q = np.zeros_like(cell_data, dtype=float)
    res = np.array([0, 0, 0, 0], dtype=float)

    # CALCULATE F AND G BASED ON Q
    dependents = np.zeros((len(x_list)-1, len(x_list[0])-1, 3, 4), dtype=float)
    dependents[:, :, 0, 0] = cell_data[:,:,1].copy()
    dependents[:, :, 0, 1] = cell_data[:,:,1]**2/cell_data[:,:,0] + (cell_data[:,:,3] - 0.5*cell_data[:,:,0]*(np.linalg.norm(cell_data[:,:,1:3], axis=2)/cell_data[:,:,0])**2)*(gamma-1)
    dependents[:, :, 0, 2] = cell_data[:,:,1]*cell_data[:,:,2]/cell_data[:,:,0]
    dependents[:, :, 0, 3] = (cell_data[:,:,3] + (cell_data[:,:,3] - 0.5*cell_data[:,:,0]*(np.linalg.norm(cell_data[:,:,1:3], axis=2)/cell_data[:,:,0])**2)*(gamma-1))*(cell_data[:,:,1]/cell_data[:,:,0])
    dependents[:, :, 1, 0] = cell_data[:,:,2].copy()
    dependents[:, :, 1, 1] = cell_data[:,:,1]*cell_data[:,:,2]/cell_data[:,:,0]
    dependents[:, :, 1, 2] = cell_data[:,:,2]**2/cell_data[:,:,0] + (cell_data[:,:,3] - 0.5*cell_data[:,:,0]*(np.linalg.norm(cell_data[:,:,1:3], axis=2)/cell_data[:,:,0])**2)*(gamma-1)
    dependents[:, :, 1, 3] = (cell_data[:,:,3] + (cell_data[:,:,3] - 0.5*cell_data[:,:,0]*(np.linalg.norm(cell_data[:,:,1:3], axis=2)/cell_data[:,:,0])**2)*(gamma-1))*(cell_data[:,:,2]/cell_data[:,:,0])
    
    # CALCULATE PRIMITIVES BASED ON Q
    primitives = np.zeros((len(x_list) - 1, len(x_list[0])-1, 3), dtype=float) # useful primitives: temperature, speeds of sound, pressure
    primitives[:,:,0] = (1/(cv))*(np.divide(cell_data[:,:,3], cell_data[:,:,0]) - (0.5)*np.pow(np.divide(np.linalg.norm(cell_data[:,:,1:3], axis=2), cell_data[:,:,0]), 2))
    primitives[:, :, 1] = np.sqrt(primitives[:, :, 0]) * np.sqrt(gamma*R) # a = sqrt(gamma*R*T)
    primitives[:, :, 2] = np.multiply(cell_data[:, :, 0], primitives[:, :, 0]) * R # p=rho*R*T

    # fs and gs, averaged
    fs = np.zeros((len(cell_data), len(cell_data[0]),    4,       4), dtype=float)
    gs = np.zeros_like(fs, dtype=float)

    # precalculate areas
    areas = np.array([[0.5*((x_list[j+1, i+1] - x_list[j, i])*(y_list[j+1, i] - y_list[j, i+1])  -  (y_list[j+1,i+1] - y_list[j,i])*(x_list[j+1,i] - x_list[j,i+1])) for i in range(len(x_list[0])-1)] for j in range(len(x_list)-1)])
    
    # __ of east, north, west, south
    # calculate x and y compoment deltas, counterclockwise (outward pointing normal)
    # cols (), rows (), delta_ys (4), delta_xs (4)
    del_y_east = lambda j,i: (y_list[j+1,i+1] - y_list[j,i+1]);       del_x_east = lambda j,i: (x_list[j+1,i+1] - x_list[j,i+1])   # delta x,y in east
    del_y_north = lambda j,i: y_list[j+1,i] - y_list[j+1,i+1];      del_x_north = lambda j,i: (x_list[j+1,i] - x_list[j+1,i+1]) # delta x,y in north
    del_y_west = lambda j,i: (y_list[j,i] - y_list[j+1,i]);           del_x_west = lambda j,i: (x_list[j,i] - x_list[j+1,i]  )         # delta x,y in west
    del_y_south = lambda j,i: y_list[j,i+1] - y_list[j,i];          del_x_south = lambda j,i: (x_list[j,i+1] - x_list[j,i])         # delta x,y in south
    cell_wall_distances = np.array([[[[del_y_east(j,i), del_y_north(j,i), del_y_west(j,i), del_y_south(j,i)],
                                   [del_x_east(j,i), del_x_north(j,i), del_x_west(j,i), del_x_south(j,i)]] 
                                   for i in range(len(x_list[0])-1)] for j in range(len(x_list)-1)])
    # negate this to make it inward pointing normal


    # lengths of cell walls
    lengths = np.sqrt(np.power(cell_wall_distances[:, :, 0, :], 2) + np.power(cell_wall_distances[:, :, 1, :], 2))
    # ainput(f'🛈 the lengths have shape {np.shape(lengths)}')

    # calculate cell wall normals
    rot90 = np.array([[0.0, -1.0], [1.0, 0.0]])
    delt_to_norm = lambda x: np.matmul(rot90, np.flip(x))/np.linalg.norm(x)
    norms = np.array([[[delt_to_norm(cell_wall_distances[j,i,:,0]),
                         delt_to_norm(cell_wall_distances[j,i,:,1]),
                           delt_to_norm(cell_wall_distances[j,i,:,2]),
                             delt_to_norm(cell_wall_distances[j,i,:,3])]
              for i in range(len(x_list[0])-1)] for j in range(len(x_list)-1)])

    # each index has one subtracted because for n nodes there are n-1 cells

    
    def get_p(j,i):
        # return ((get_rho(j,i)) * get_T(j,i)) / (T_ref) # calculated directly from other nondimensional values
        j = ind(j)
        i = ind(i)
        return max(0.0, primitives[j,i,2]) # calculated by getting dimensional pressure (rho_d * R * T_d) and nondimensionalizeing as recommended by cizmas (p / rho u_ref**2)
    
    # dissipation calculators
    def ind(k):
        """indexify an index"""
        assert np.isclose(k%1, 0), f"Index issue - irl your indeces can only be integers, not {k}"
        return int(np.round(k))

    def lamb(j, i, og=-1/2):
        """ becomes lamb_xi or lamb_eta depending on what axis the off-integer value is in \n
            off integer in i  -  lamb xi  -  north velocity
            off integer in j  -  lamb eta -  east velocity

            og -1/2 means the cell to the left of your wall
            og +1/2 means the cell to the right of your wall
        input: off-integer in EITHER j or i \n
        returns: the absolute wall-normal velocity component, with respect to the wall corresponding to the half-index """

        mults =  np.isclose(np.array([j%1, i%1]), np.array([0.5, 0.5]), rtol=0.05).astype(np.float64)

        # og determines whether to get velocity from left or right cell, negative meaning left cell, positive meaning right cell

        # bottom/left node (aka, truncated)
        j_base = math.floor(j) # 5
        i_base = math.floor(i) # 4

        # whether to round up or not - always round up the one with the half index
        j_round = 1 if abs(j - j_base - 0.5)<1e-3 else 0
        i_round = 1 if abs(i - i_base - 0.5)<1e-3 else 0

        # cell index
        #    - move up or down depending on og (minus moves you back)
        #    - move i or j depending on who is off-integer
        jc = ind(j + og*j_round )
        ic = ind(i + og*i_round)
        V = cell_data[jc, ic, 1:3] / cell_data[jc, ic, 0] # vector velocity for right/top cell

        # sign of the normal does not matter, only the direction
        # what if the sign of the normal DOES matter??
        mults =  np.isclose(np.array([j%1, i%1]), np.array([0.5, 0.5]), rtol=0.05).astype(np.float64)
        n = norms[j_base, i_base, np.argmax(mults), :2] # wall's unit normal

        a_cell = primitives[jc, ic, 1].copy()

        vn = np.dot(n, V)

        # if i>20 and j >= 1:
        #     print(f'EIGENVALUE requested at {j,i}, {og}')
        #     print(f'measured between the nodes (j,i = {jc, ic} and (j,i = {jc, ic})')
        #     print(f'delta x, delta_y = {delta_x, delta_y}')
        #     print(f'mults: {mults}')
        #     print(f'norms for {jc, ic}: {n}')
        #     print(f'{['east', 'north'][np.argmax(mults)]}')
        #     print(f'normal: {n} with angle of {np.rad2deg(np.atan2(n[1], n[0]))}')
        #     print(f'measured vevlocity at {jc, ic}')
        #     print(f'vn_abs: {vn}\n\
        #         a: {a_cell}')
        #     print(f'lambda: {vn+a_cell}')
        #     plt.quiver(0.5*(x_list[j_base, i_base]+x_list[j_base+j_round, i_base+i_round]), 0.5*(y_list[j_base, i_base]+y_list[j_base+j_round, i_base+i_round]), n[0], n[1])
        #     plt.scatter(x_list[jc,ic]+0.5*delta_x, y_list[jc,ic]+0.5*delta_y)
        #     mesh_plot(x_list, y_list, 'u')

        return vn + a_cell

    def switch2_eta(jc,i, j):
        adder = 1 if jc>j else -1

        j = ind(j); i = ind(i)

        num = abs(primitives[j+1+adder,i,2] - 2*primitives[j+adder,i,2] + primitives[j-1+adder,i,2])
        den =     primitives[j+adder,i+1,2] + 2*primitives[j+adder,i,2] + primitives[j+adder,i-1,2] + 1e-3
        res1 = num/den

        num = abs(primitives[j+1,i,2] - 2*primitives[j,i,2] + primitives[j-1,i,2])
        den =     primitives[j,i+1,2] + 2*primitives[j,i,2] + primitives[j,i-1,2] +  1e-3
        res2 = num/den

        return nu_2*0.5*(res1+res2)

    def switch2_xi(j,ic, i):
        # j is on integer, and ic is off-integer. i is the home cell.
        # returns the difference centered about the home cell averaged with either the right or the left depending on if it was +1/2 or -1/2
        
        # the switch at the face is the average of the switch at the cells - cizmas
        adder = 1 if ic>i else -1
        j = ind(j); i = ind(i)
 

        num = abs(primitives[j,i+1+adder,2] - 2*primitives[j,i+adder,2] + primitives[j,i-1+adder,2])
        den =     primitives[j,i+1+adder,2] + 2*primitives[j,i+adder,2] + primitives[j,i-1+adder,2] + 1e-3
        res1 = num/den

        num = abs(primitives[j,i+1,2] - 2*primitives[j,i,2] + primitives[j,i-1,2])
        den =     primitives[j,i+1,2] + 2*primitives[j,i,2] + primitives[j,i-1,2] + 1e-3
        res2 = num/den

        return nu_2*0.5*(res1+res2)


    def switch4_xi(jc,ic, i):
        #jc is on integer
        return np.max([0.0, (nu_4)-switch2_xi(jc, ic, i)])
    
    def switch4_eta(jc,ic, j):
        return np.max([0.0, (nu_4)-switch2_eta(jc, ic, j)])

    def enforce_BCs():
        # farfield (inlet boundary condition)
        for j in range(len(x_list)-1): # for all the rows, -1 because n-1 cells for n nodes
            for i in [0, 1]: # for the two ghost cells going in (same inlet condition)
                cell_data[j, i, :] =  np.array([rho_infty,
                                            mach*a*rho_infty,
                                            0.0,
                                            rho_infty*E_inlet], dtype=float)
            
            # ZERO gradient at outlet EXCEPT FOR ENERGY
            # calculate state from zero gradient
            # q_bc
            cell_data[j, -2:, :3] = 2*cell_data[j,-3, :3] - cell_data[j,-4, :3]

            # energy
            V = cell_data[j, -2, 1:3]/cell_data[j, -2, 0]
            V_mag = np.sqrt(np.dot(V, V))

            # exhaust BC - same pressure as air coming in            
            cell_data[j,len(cell_data[0])-2:, 3] = p_infty / (gamma - 1) + 0.5*cell_data[j, -2, 0]*V_mag**2 # = rho_E_out.copy()
            # get the normal to the border wall

        for m in [-1, 1]: # for top, bottom
            for Aa in m*np.array([[0, -1], [1, -2]]): # (wall-border cell, wall-border ghost), (non-wall-border cell, non-wall-border ghost)
                border_j = 5*int(m>=0)-3 # -3 for the top, 2 for the bottom (nodes)
                domain_j, ghost_j = Aa[:] + border_j # cells
                for i in range(np.shape(x_list)[1]-1): # iterate down axis 1, -1 because index 55 is out of bounds for axis 1 with size 55

                    # get the normal to the border wall
                    delta_x = x_list[border_j, i+1] - x_list[border_j, i]
                    delta_y = y_list[border_j, i+1] - y_list[border_j, i]
                    norm_factor = 1/ (delta_x**2 + delta_y**2)**0.5 # to normalize
                    t_wall = np.array([delta_x, delta_y]) * norm_factor
                    n_wall = m*np.array([-delta_y, delta_x]) * norm_factor # multiply by n to make sure the normals point inward

                    assert np.isclose(np.linalg.norm(t_wall), 1), np.linalg.norm(t_wall)
                    assert np.isclose(np.linalg.norm(n_wall), 1), np.linalg.norm(n_wall)
                    assert np.isclose(np.dot(n_wall, t_wall), 0)

                    # get the momentum for the referenced cell, and reflect the normal component and preserve the tangent compoonent
                    reference_V = cell_data[domain_j, i, 1:3] / cell_data[domain_j, i, 0]
                    reflected_V = t_wall * (np.dot(reference_V, t_wall)) - n_wall * (np.dot(reference_V, n_wall))
                    # reflected_P = reference_P @ rotator(-2*angle)
                    cell_data[ghost_j, i, :] = np.array([cell_data[domain_j, i, 0],
                                                        reflected_V[0]*cell_data[domain_j, i, 0],
                                                            reflected_V[1]*cell_data[domain_j, i, 0],
                                                            cell_data[domain_j, i, 3]])
                    assert np.isclose(np.linalg.norm(reference_V), np.linalg.norm(reflected_V))
                    # print(f'angle of flow with wall: {angle}')
                    net_V = cell_data[ghost_j, i, 1:3] + cell_data[domain_j, i, 1:3] # ensures boundary conditions are flawless
                    assert np.isclose(np.arctan2(net_V[1], net_V[0]), np.arctan2(t_wall[1], t_wall[0])), f'failed to implement BCs: momentum: {net_V}, wall tangent: {t_wall}'
                    # also throws an error when there is reversed flow
                    # print(f'domain j, ghost j, border j  = {domain_j, ghost_j, border_j}')
                    if not np.isclose(np.linalg.norm(cell_data[domain_j, i, 1:3]), np.linalg.norm(cell_data[ghost_j, i, 1:3])):
                        print(f'mags supposed to be equal: {np.linalg.norm(cell_data[domain_j, i, 1:3]), np.linalg.norm(cell_data[ghost_j, i, 1:3])}')
                        print(f'{cell_data[domain_j, i, 1:3]} - {cell_data[ghost_j, i, 1:3]}')
                    assert np.isclose(np.linalg.norm(cell_data[domain_j, i, 1:3]), np.linalg.norm(cell_data[ghost_j, i, 1:3]))

                   
                    # *special* energy calculation - possible bugs?
                    new_rho_E = (get_p(domain_j, i) / (gamma - 1)) + 0.5 * np.dot(cell_data[ghost_j, i, 1:3], cell_data[ghost_j, i, 1:3]) / cell_data[ghost_j, i, 0]
                     # assert np.isclose(new_rho_E, cell_data[ghost_j, i, 3] ), F"{new_rho_E}!={cell_data[ghost_j, i, 3]}"
                    cell_data[ghost_j, i, 3]  = float(new_rho_E)

    # ========================= NOTES =========================
    # NOTE: when iterating over cells instead of nodes:
    #   the bottom and left nodes are going to be i,j
    #   the top and right nodes are going to be i+1,j+1
    #   make sure to iterate up until the n-1th term for a list with n entries

    # NOTE: values at the borders are calculated with averages

    # ========================= BOUNDARY CONDITIONS =========================

    # set up the boundary conditions here because you want to re-enforce them at every time step

    # NOTE: similarly to the mesher, I want to have them outside the main loop so that I can test the funcitonaltiy independent of the main solver
    # later, I can integrate them into the main ij loop for better speeds

    # check cizmas p.167
    enforce_BCs()
    while True: 
        ag_res = 0.0 # average residual for the current iteration
        it_number += 1 # my iteration number just for my own bookkeeping
        print(it_number, end='      ')

        enforce_BCs()

        # calculate some useful primitives: temperature, speed of sound, pressure
        primitives[:,:,0] = (1/(cv))*(np.divide(cell_data[:,:,3], cell_data[:,:,0]) - (0.5)*np.pow(np.divide(np.linalg.norm(cell_data[:,:,1:3], axis=2), cell_data[:,:,0]), 2))
        primitives[:, :, 1] = np.sqrt(primitives[:, :, 0]) * np.sqrt(gamma*R) # a = sqrt(gamma*R*T)
        primitives[:, :, 2] = np.multiply(cell_data[:, :, 0], primitives[:, :, 0]) * R # p=rho*R*T

        enforce_BCs() # update boundary conditions after its dependencies get updated

        # update all Fs and Gs after boundary conditions have been enforced
        dependents[:, :, 0, 0] = cell_data[:,:,1].copy()
        dependents[:, :, 0, 1] = cell_data[:,:,1]**2/cell_data[:,:,0] + (cell_data[:,:,3] - 0.5*cell_data[:,:,0]*(np.linalg.norm(cell_data[:,:,1:3], axis=2)/cell_data[:,:,0])**2)*(gamma-1)
        dependents[:, :, 0, 2] = cell_data[:,:,1]*cell_data[:,:,2]/cell_data[:,:,0]
        dependents[:, :, 0, 3] = (cell_data[:,:,3] + (cell_data[:,:,3] - 0.5*cell_data[:,:,0]*(np.linalg.norm(cell_data[:,:,1:3], axis=2)/cell_data[:,:,0])**2)*(gamma-1))*(cell_data[:,:,1]/cell_data[:,:,0])
        
        dependents[:, :, 1, 0] = cell_data[:,:,2].copy()
        dependents[:, :, 1, 1] = cell_data[:,:,1]*cell_data[:,:,2]/cell_data[:,:,0]
        dependents[:, :, 1, 2] = cell_data[:,:,2]**2/cell_data[:,:,0] + (cell_data[:,:,3] - 0.5*cell_data[:,:,0]*(np.linalg.norm(cell_data[:,:,1:3], axis=2)/cell_data[:,:,0])**2)*(gamma-1)
        dependents[:, :, 1, 3] = (cell_data[:,:,3] + (cell_data[:,:,3] - 0.5*cell_data[:,:,0]*(np.linalg.norm(cell_data[:,:,1:3], axis=2)/cell_data[:,:,0])**2)*(gamma-1))*(cell_data[:,:,2]/cell_data[:,:,0])
        
        # calculate fs and gs at the cells

        # roll -1 in axis 1 to see your neighbor to the east (right)
        # roll 1 in axis 0 to see your neighbor to the north,
        # roll 1 in axis 1 to see you neigbber to the west (left)
        # roll -1 in axis 0 to see your neighbor to the south
        rdirec = [-1, 1, 1, -1]  
        axisnum = [1, 0, 1, 0]
        f_direc = [0, 1, 2, 3]
        for rd, ax, he in zip(rdirec, axisnum, f_direc):
            # you must roll -1 to see your neighbor to the north, and +1 to see your neightbor in the south, etc
            fs[:, :, he, :] = (np.roll(dependents[:, :, 0, :], rd, axis=ax) + dependents[:, :, 0, :])*0.5
            gs[:, :, he, :] = (np.roll(dependents[:, :, 1, :], rd, axis=ax) + dependents[:, :, 1, :])*0.5


        # =============== INNNER EQ =========================
        for j in range(2, len(x_list)-2 - 1):
            for i in range(2, len(x_list[0])-2 - 1):
                 
                assert np.shape(fs[j,i, :, :]) == (4, 4), f'shape of fs is {np.shape(fs[j,i, :, :])}'
                assert np.shape(cell_wall_distances[j,i,0,:]) == (4,), f'shape of cellwalls is {np.shape(cell_wall_distances[j,i,0,:])}'
                
                res = np.subtract(fs[j,i, :, :].T @ cell_wall_distances[j,i,0,:],
                            gs[j,i, :, :].T @ cell_wall_distances[j,i,1,:]).copy()

                # use JST dissipation term to reduce wiggles - https://i.imgur.com/pdUVsjX.png
                if (it_number)%4 == 0:
 
                    term1 = switch2_xi(j,i+1/2,i+1) * lamb(j,i+1/2,og=-1/2) * lengths[j,i,0] * (cell_data[j,i+1] - cell_data[j,i]) - \
                            switch2_xi(j,i-1/2,i-1) * lamb(j,i-1/2,og=+1/2) * lengths[j,i,2] * (cell_data[j,i] - cell_data[j,i-1])
                    
                    term2 = switch2_eta(j+1/2,i,j+1) * lamb(j+1/2,i,og=-1/2) * lengths[j,i,1] * (cell_data[j+1,i] - cell_data[j,i]) - \
                            switch2_eta(j-1/2,i,j-1) * lamb(j-1/2,i,og=+1/2) * lengths[j,i,3] * (cell_data[j,i] - cell_data[j-1, i])
                    
                    term3 = switch4_xi(j,i+1/2,i+1) * lamb(j,i+1/2,og=-1/2) * lengths[j,i,0] * (cell_data[j,i+2] - 3*cell_data[j,i+1] + 3*cell_data[j,i] - cell_data[j,i-1]) - \
                            switch4_xi(j,i-1/2,i-1) * lamb(j,i-1/2,og=+1/2) * lengths[j,i,2] * (cell_data[j,i+1] - 3*cell_data[j,i] + 3*cell_data[j,i-1] - cell_data[j,i-2])

                    term4 = switch4_eta(j+1/2,i,j+1) * lamb(j+1/2,i,og=-1/2) * lengths[j,i,1] * (cell_data[j+2,i] - 3*cell_data[j+1, i] + 3*cell_data[j,i] - cell_data[j-1,i]) - \
                            switch4_eta(j-1/2,i,j-1) * lamb(j-1/2,i,og=+1/2) * lengths[j,i,3] * (cell_data[j+1,i] - 3*cell_data[j,i] + 3*cell_data[j-1,i] - cell_data[j-2,i])
                    
                    dependents[j,i,2,:] = term1 + term2 - term3 - term4

                # RK4 integration parameters - entry by entry multiplication of eigs and lens, then sum them all
                delta_t = 2 * CFL * areas[j,i] / np.sum(np.multiply(np.array([lamb(jt,it, og=((j-jt)+(i-it))) for jt,it in zip([j, j+1/2, j, j-1/2], [i+1/2, i, i-1/2, i])]), lengths[j,i]))

                # calculate new state
                new_q[j,i,:] = cell_data[j,i,:] - (alphas[(it_number-1)%4]*delta_t / areas[j,i]) * (res - dependents[j,i,2,:]) # RK iteration with CFL

                for k in range(4):
                    if j == 21 and i == 2:
                        print(f'new_q[{k}] = {cell_data[j,i,:]} - {(alphas[(it_number-1)%4]*delta_t / areas[j,i])} * ({res} - {dependents[j,i,2,:]})')
                # laplacian dissipation
                # new_q[j,i,:] += 1e-3*(cell_data[j, i+1, :] + cell_data[j, i-1, :] + cell_data[j+1, i, :] + cell_data[j-1, i, :] - 4*cell_data[j,i, :])
                # input(f'residual at ({j,i})={res}')
                # residual calculations
                ag_res += abs(np.linalg.norm(res))

        # do this regardless
        cell_data[2:-2,2:-2,:] = new_q[2:-2, 2:-2, :].copy()
        enforce_BCs()
        
        # if it_number%300==0:
        #     total_res -= cell_data
        #     plot_scalars(x_list, y_list, total_res[:, :, 3], caption='total_res in energy', verbose=True)

        # status updates
        # prevent nan residuals
        assert not np.isnan(ag_res), "residual cannot be nan"

        ag_res_list.append(ag_res/(len(x_list)*len(x_list[0])))
        print(f'Residual Magnitude={ag_res_list[-1]:.3E}\t\t System v: {np.sum(cell_data[:, :, 2]/cell_data[:, :, 0]):.2E}\t\t u: {np.sum(cell_data[:, :, 1]/cell_data[:, :, 0]):.6f}\t\t rho: {np.sum(cell_data[:, :, 0]):.6f}\t\t E: {np.sum(cell_data[:, :, 3]/cell_data[:, :, 0]):.6f}')
        tot_mass_list.append(np.dot(areas.ravel(), (cell_data[:, :, 0]).ravel()))
        avg_u_list.append(np.dot(areas.ravel(), (cell_data[:, :, 1]/cell_data[:, :, 0]).ravel()))
        avg_v_list.append(np.dot(areas.ravel(), (cell_data[:, :, 2]/cell_data[:, :, 0]).ravel()))
        tot_energy_list.append(np.dot(areas.ravel(), (cell_data[:, :, 3]/cell_data[:, :, 0]).ravel()))


        # major status updates (every 300)
        if int(it_number%300) == 0:

            plt.plot(np.linspace(0, len(ag_res_list), num=len(ag_res_list)), ag_res_list)
            plt.yscale('log');  plt.title(f'Residual Magnitude')
            plt.savefig(f'{savedir}current_residuals - jmax,imax={len(x_list)}x{len(x_list[0])}, M={mach}.png'); plt.clf()

            # debugging - net mass, energy, density, etc
            plt.plot(np.linspace(0, len(tot_mass_list), num=len(avg_u_list)), avg_u_list)
            plt.yscale('linear');  plt.title(f'System Mass (kg) as a function of iterations')
            plt.savefig(f'{savedir}system mass - jmax,imax={len(x_list)}x{len(x_list[0])}, M={mach}.png'); plt.clf()

            plt.plot(np.linspace(0, len(avg_u_list), num=len(avg_u_list)), avg_u_list)
            plt.yscale('linear');   plt.title(f'system bulk u')
            plt.savefig(f'{savedir}system bulk u - jmax,imax={len(x_list)}x{len(x_list[0])}, M={mach}.png'); plt.clf()

            plt.plot(np.linspace(0, len(avg_v_list), num=len(avg_v_list)), avg_v_list)
            plt.yscale('linear');   plt.title(f'system bulk v')
            plt.savefig(f'{savedir}system bulk v- jmax,imax={len(x_list)}x{len(x_list[0])}, M={mach}.png'); plt.clf()

            plt.plot(np.linspace(0, len(tot_energy_list), num=len(tot_energy_list)), tot_energy_list)
            plt.yscale('linear');   plt.title(f'system total energy')
            plt.savefig(f'{savedir}system total energy - jmax,imax={len(x_list)}x{len(x_list[0])}, M={mach}.png'); plt.clf()


            plot_scalars(x_list, y_list, cell_data[:, :, 0], caption=f'density - jmax,imax={len(x_list)}x{len(x_list[0])}, M={mach}') # plot densities (debugging)
            plot_scalars(x_list, y_list, cell_data[:, :, 2], caption=f'F2 - jmax,imax={len(x_list)}x{len(x_list[0])}, M={mach}') # plot densities (debugging)

            plot_scalars(x_list, y_list, np.divide(cell_data[:, :, 3], cell_data[:, :, 0]), caption=f'energy - jmax,imax={len(x_list)}x{len(x_list[0])}, M={mach}') # plot densities (debugging)
            plot_scalars(x_list, y_list, dependents[:, :, 2, 3], caption=f'dissipations - energy jmax,imax={len(x_list)}x{len(x_list[0])}, M={mach}') # plot densities (debugging)

            # x velocities and y velocities
            us = np.divide(cell_data[:, :, 1].ravel(), cell_data[:, :, 0].ravel()) # get u velocities
            vs = np.divide(cell_data[:, :, 2].ravel(), cell_data[:, :, 0].ravel()) # get v velocities
            machs = np.power(np.power(us, 2) + np.power(vs, 2), 0.5) / a
            plot_scalars(x_list, y_list, machs, caption=f'Machs - jmax,imax={len(x_list)}x{len(x_list[0])}, M={mach}')
            #  plot_quivers(x_list,y_list, us, vs, caption=f'velocities', verbose=True)

            ravelled_size = np.shape(cell_data)[0] * np.shape(cell_data)[1] * np.shape(cell_data)[2]
            longer_size = int(max(ravelled_size, len(ag_res_list)))
            
            # create a DataFrame
            data_padded = cell_data.ravel().tolist() + [np.nan] * (longer_size - ravelled_size)
            residuals_padded = ag_res_list + [np.nan] * (longer_size - len(ag_res_list))
            df = pd.DataFrame({'data': data_padded,'convergence': residuals_padded})
            # save work no matter what
            try:
                df.to_csv(f'{savedir}j_max={len(x_list)-1},i_max={len(x_list[0])-1},M={round(mach, 1)}.csv')
            except KeyboardInterrupt:
                print(f"I'm writing let me finish.")
                df.to_csv(f'{savedir}j_max={len(x_list)-1},i_max={len(x_list[0])-1},M={round(mach, 1)}.csv')
                exit()
            del df # clean up so there are no memory overflows

           # print out the force on the bottom bump as well as a mach plot
            corner_pts = [find_arg_nearest(x_list[0], 2.0), find_arg_nearest(x_list[0], 3.0)]
            Fx = 0.0
            Fy = 0.0
            for i in range(corner_pts[0], corner_pts[1]):
                j = 2
                delta_y = (y_list[ind(j+1), ind(i+1)] - y_list[ind(j+int(np.round(j-j+epsilon))), ind(i+int(np.round(i-i+epsilon)))])
                delta_x = (x_list[ind(j+1), ind(i+1)] - x_list[ind(j+int(np.round(j-j+epsilon))), ind(i+int(np.round(i-i+epsilon)))])
                n = np.array([-delta_y, delta_x]) / (delta_y**2 + delta_x**2)**0.5
                F = get_p(j,i) * lengths[2, i, 3] * n
                Fx += F[0]
                Fy += F[1]
            msg = f'X force: {Fx} N\n\
                    Y force: {Fy} N'
            # print(msg)
            with open(f'{savedir}magic number- M={mach}, ji={len(x_list),len(x_list[0])}.txt', 'w') as my_glorious_result:
                my_glorious_result.write(msg)

            if it_number >=5000 and False:
                exit()
                if np.average(ag_res_list[-20:-1]) <= 0.01:
                    print(f'Converged to within 1%')
                    break

    return 0


# main sequence
def main():
    # ensure directory exists
    if not os.path.exists(f'{savedir}'):
        os.makedirs(f'{savedir}')


    machs = [0.3]#, 0.3, 0.7]
    for M in machs:
        print(f'== MACH {M} ==')
        # shapes = [(11,51)]#, (21,101), (41,201)]
        # shapes = [(11,51)]#, (21,101), (11, 51)]
        shapes = [(10, 50)]

        for shape in shapes:
            # get the mesh
            try:
                df = pd.read_csv('saved_meshes.csv')
            except FileNotFoundError:
                df = mesh(shape)
                df.to_csv('saved_meshes.csv')
            except pd.errors.EmptyDataError:
                df = mesh(shape)
                df.to_csv('saved_meshes.csv')
            try:
                x_list = df[f'{shape}-x']
                y_list = df[f'{shape}-y']
                print(f'found pre-calculated mesh')
            except KeyError:
                df_addition = mesh(shape)
                df = pd.concat([df, df_addition], axis=1)
                df.to_csv('saved_meshes.csv')
            
            x_list = trim_dataframe(df[f'{shape}-x']).to_numpy().reshape(shape)
            y_list = trim_dataframe(df[f'{shape}-y']).to_numpy().reshape(shape)
            # print(f"{shape} mesh has average skewness of: {mesh_plot(x_list, y_list, 'FINAL MESH', return_squishedness=True)}")
            
            # run the simulation and generate results
            run(x_list, y_list, M)

    print(f'Process finished with exit code 0')
    exit()

# run the script
if __name__ == "__main__":
    main()