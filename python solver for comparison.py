import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pyvista as pv
from numpy import float64
import os
epsilon = 1e-5

class Mesh:
    def __init__(self, dimensions):
        """Input the dimensions of the physical mesh (number of core cells tall by number of core cells wide)"""
        self.core_dimensions = dimensions
    def generate(self):
        """Generates a mesh in un-ghostified form."""
        min_residual=1e-6

        # Resolution:
            # - resolution of under 4 causes everything to completely break
            # - resolution is nodes per unit
            # - change the resolution to prove mesh-independence of solution

        # min residual:
            # - the average amount that is changed by the iterator
            # - default is 1e-6

        def find_arg_nearest(array, value):
            """inputs: the array you want to search, value that you want to seach for
                outputs: the index of the value closest to the value you searched for"""
            idx = np.searchsorted(array, value, side="left")
            if idx > 0 and (idx == len(array) or np.fabs(value - array[idx-1]) < np.fabs(value - array[idx])):
                return idx-1
            else:
                return idx
            
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

        # define number of nodes
        # +1 because number of cells-->number of nodes
        i_max = self.core_dimensions[0]+1 # y direction (axis 0)
        j_max = self.core_dimensions[1]+1 # x direction (axis 1)

        # generate xi-eta grid
        # xi (x)
        xi_list = np.zeros((i_max, j_max), dtype=float) # create J_max rows of I_max entries each
        for i in range(i_max): # for every row in the J_max steps
            for j in range(j_max): # generate a row with I_max steps
                xi_list[i,j] = (i) / (i_max-1)

        # eta (y)
        eta_list = np.zeros((i_max, j_max), dtype=float) # create __ rows of ___ entries each
        for i in range(i_max): # for every row,
            for j in range(j_max): # generae a row.
                eta_list[i,j] = (j) / (j_max-1)
        
        # xi_list, eta_list = np.meshgrid(np.linspace(0, 1, num=j_max), np.linspace(0,1, num=i_max))

        # test plot
        plt.title(f'Algebraic Grid in ξ-η Space'); plt.plot(eta_list, xi_list, marker='.', color='k', linestyle='none'); plt.xlabel('eta, or j'); plt.ylabel('xi, or i'); plt.show()
        
        # algebraic grid
        x_list = np.zeros((i_max, j_max), dtype=float)
        y_list = x_list.copy() # pre-allocate memory. Remember - zeroes works by generating __ rows of ___ length
        for i in range(i_max):
            for j in range(j_max):
                try:
                    x_list[i, j] = x_l + eta_list[i, j] * (x_r - x_l) # linearly scale x up
                except:
                    print(f'bad x, {[i, j]} out of bounds for either x list: {np.shape(x_list)} or xi_list: {np.shape(xi_list)}')
                    print(x_list)
                    exit()
                try:
                    x = x_list[0, j] # COULD take advantage of the fact that the domain on the x scales from 0 to 1 and just say xi, but eh
                    y_list[i, j] = y_lower(x) + xi_list[i,j] * (y_upper(x) - y_lower(x))
                except:
                    print(f'bad y')
                    print(y_list)
                    exit()

        plt.title(f'Algebraic Grid in x-y Space'); plt.plot(x_list, y_list, marker='.', color='k', linestyle='none'); plt.show()

        # preliminary calcs
        delta_xi = 1/(j_max)
        delta_eta = 1/(i_max) # this is because xi and eta are defined to be on the interval {0<x<1}. In python if your list has max index of 3, then you have 2 divisions

        # break out of the loop either when residuals hit the floor or iterations hit the ceiling
        max_iterations = 100*self.core_dimensions[0] # iteration ceiling
        # min_residual = 1e-6 # residual floor

        aggregate_residuals = np.zeros((max_iterations+1))
        print(f'Beginning a maximum of {max_iterations} iterations.\nCurrent iteration:  ', end='')

        iterations = 0
        while iterations < max_iterations:
            iterations += 1

            # fancy schmancy status
            print('\b'*len(str(iterations-1)), end='', flush=True)
            print(iterations, end='', flush=True)

            corner_pts = [find_arg_nearest(x_list[0, :], 2.0), find_arg_nearest(x_list[0, :], 3.0)] # define the corner points and the bool function (for simplicity)
            def corner_point(i,j):
                if i in [0, i_max-1] and j in corner_pts:
                    return True
                else:
                    return False

            # sets the iterations for i,j
            iter_set = np.array([[(i, j) for j in range(j_max)] for i in range(i_max)]).reshape(i_max*j_max, 2)

            iter_residual = 0
            for i,j in iter_set:
                # corner points
                if (corner_point(i, j)):
                    if i == 0 and j < j_max/2:
                        new_x_val = 2.0
                        new_y_val = 0.0
                        y_list[i, j] = new_y_val
                    elif i == i_max - 1 and j < j_max/2:
                        new_x_val = 2.0
                        new_y_val = 1.0
                        y_list[i, j] = new_y_val
                    elif i == 0 and j > j_max/2:
                        new_x_val = 3.0
                        new_y_val = 0.0
                    elif i == i_max-1 and j > j_max/2:
                        new_x_val = 3.0
                        new_y_val = 1.0

                # boundary conditions
                elif i in [0, i_max - 1]: # for top and bottom
                    if x_list[i, j] < 2 or x_list[i, j] > 3:
                        try: 
                            new_x_val = x_list[i+1, j]
                        except:
                            new_x_val = x_list[i-1, j]
                        new_x_val = x_list[i,j] # override for corner points
                    if x_list[i, j] < 3 and x_list[i, j] > 2:
                        try:
                            new_x_val = x_list[i+2, j] + ((y_lower(x_list[i, j-1]) - y_lower(x_list[i, j+1]))/(x_list[i,j-1] - x_list[i, j+1])) * delta_xi * x_r 
                            new_y_val = y_lower(new_x_val)
                            iter_residual += abs(y_list[i,j] - new_y_val)
                            y_list[i,j] = new_y_val
                        except:
                            new_x_val = x_list[i-2,j] + ((y_lower(x_list[i, j-1]) - y_lower(x_list[i,j+1]))/(x_list[i,j-1] - x_list[i,j+1])) * delta_xi * x_r 
                            new_y_val = y_upper(new_x_val)
                            iter_residual += abs(y_list[i,j] - new_y_val)
                            y_list[i,j] = new_y_val
                    iter_residual += abs(x_list[i, j] - new_x_val)
                    x_list[i,j] = new_x_val
                elif j in [0, j_max - 1]: # for left and right
                    if j == 0:
                        new_y_val = y_list[i, j+1]
                        iter_residual += abs(y_list[i, j] - new_y_val)
                    elif j == j_max - 1:
                        new_y_val = y_list[i, j - 1]
                    iter_residual += abs(y_list[i, j] - new_y_val)
                    y_list[i,j] = new_y_val
                
                else: # general inside equations
                    alpha = 1/(4*delta_eta**2) * ((x_list[i,j+1] - x_list[i,j-1])**2 + (y_list[i,j+1] - y_list[i,j-1])**2) # see cizmas p.33
                    beta = -1/(4*delta_eta*delta_xi) * ((x_list[i+1, j] - x_list[i-1, j])*(x_list[i, j+1] - x_list[i, j-1])   +   (y_list[i+1, j] - y_list[i-1, j])*(y_list[i,j+1] - y_list[i,j-1]))
                    gamma = 1/(4*delta_xi**2) * ((x_list[i+1,j] - x_list[i-1,j])**2 + (y_list[i+1,j] - y_list[i-1,j])**2)

                    k1 = alpha / delta_xi**2
                    k2 = -beta / (2*delta_xi*delta_eta)
                    k3 = gamma / delta_eta**2

                    new_x_val = (0.5/(k1+k3)) * (k2*x_list[i-1, j-1] + k3*x_list[i, j-1] -k2*x_list[i+1, j-1]+ k1*x_list[i-1, j]+ k1*x_list[i+1, j] -k2*x_list[i-1, j+1] + k3*x_list[i, j+1] + k2*x_list[i+1, j+1])
                    x_list[i,j] = new_x_val
                    iter_residual += abs(x_list[i, j] - new_x_val)

                    new_y_val = (0.5/(k1+k3)) * (k2*y_list[i-1, j-1] + k3*y_list[i, j-1] -k2*y_list[i+1, j-1]+ k1*y_list[i-1, j]+ k1*y_list[i+1, j] -k2*y_list[i-1, j+1] + k3*y_list[i, j+1] + k2*y_list[i+1, j+1])
                    iter_residual += abs(y_list[i, j] - new_y_val)
                    y_list[i,j] = new_y_val

            aggregate_residuals[iterations] = iter_residual/(i_max * j_max)
            if aggregate_residuals[iterations] < min_residual:
                print(f'\nResidual<{min_residual} - convergence detected!')
                break
        else:
            print(f'\niteration ceiling reached')

        self.x_list = x_list
        self.y_list = y_list
    def ghostify(self):
        """appends ghost cells"""
        delta_x = np.average(np.diff(self.x_list[0, :]))
        delta_y = np.average(np.diff(self.y_list[:, 0].ravel()))

        # append in the top and bottom
        self.x_list = np.concatenate(([self.x_list[0, :]], [self.x_list[0, :]], self.x_list[:, :], [self.x_list[-1, :]], [self.x_list[-1, :]]), axis=0)
        self.y_list = np.concatenate(([self.y_list[0, :] - 2*delta_y], [self.y_list[0, :] - delta_y], self.y_list[:, :], [self.y_list[-1, :] + delta_y], [self.y_list[-1, :] + 2*delta_y]), axis=0)

        # append in the left and the right
        self.x_list = np.concatenate((self.x_list[:, 0:1] - 2*delta_x, self.x_list[:, 0:1] - delta_x, self.x_list[:], self.x_list[:,  -1:] + delta_x, self.x_list[:,  -1:] + 2*delta_x), axis=1)
        self.y_list = np.concatenate((self.y_list[:, 0:1], self.y_list[:, 0:1], self.y_list[:], self.y_list[:,  -1:], self.y_list[:,  -1:]), axis=1)
    def unghostify(self):
        """Removes ghost cells
        untested, but very simple and should work"""
        self.x_list = self.x_list[2:-2]
        self.y_list = self.y_list[2:-2]
    def save(self, fname):
        """save the mesh. For less headaches, the saved meshes must always have the ghost cells attached."""
        df = pd.DataFrame({f'x':self.x_list.ravel(),
                            f'y':self.y_list.ravel()})
        df.to_csv(fname)
        print(f'New Mesh saved to {os.path.join(os.getcwd(), fname)}')

    def read(self, fname):
        ghostified_dimensions = [i+5 for i in self.core_dimensions]
        df = pd.read_csv(fname)
        self.x_list = df[f'x'].to_numpy().reshape(ghostified_dimensions)
        self.y_list = df[f'y'].to_numpy().reshape(ghostified_dimensions)
        print(f'Mesh read from {os.path.join(os.getcwd(), fname)}')
    def plot(self):
        """Pass data in as x_list and y_list as a list of lists in the shape of the mesh"""
        self.squishedness = mesh_plot(self.x_list, self.y_list)
        plt.show() 



def mesh_plot(x_list, y_list, title='Mesh'):
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
    plt.axis('equal')
    return squishedness/(i_ubound*j_ubound)



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
    return 0
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
    # scalar_data = np.random.rand(grid.n_cells) 
    grid.cell_data[caption] = scalar_data.ravel()
    # Visualize the mesh
    
    plotter = pv.Plotter(off_screen=(not verbose))
    # plotter.add_mesh(grid, color='white', show_edges=True)
    plotter.add_mesh(grid, scalars=caption, cmap='viridis', show_edges=True)
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


def run(x_list, y_list, mach):

    loud = False    # whether to print out debugging messages

    print(f'{np.shape(x_list)}')
    
    # pre-allocate tools
    D = np.zeros(4, dtype=float) # pre-allocated dissipation
    new_q = np.zeros(4, dtype=float)

    # # f, g, and h are fluxes through the x, y, and z faces
    f = np.array([[1, 1, 1, 1],# east x flux 
                [1, 1, 1, 1], # north x flux
                [1, 1, 1, 1], # west x flux
                [1, 1, 1, 1]], dtype=float) # south x flux ()
    
    g = np.array([[1, 1, 1, 1],
                [1, 1, 1, 1], 
                [1, 1, 1, 1], 
                [1, 1, 1, 1]], dtype=float) # central Y flux - the "flux" is caluclated at the cell, and you get the actual flux values by averaging between two cells

    
    # mesh_plot(x_list, y_list, 'Ghost cells added') # don't need to see this every time LOL

    # initial gas properties
    gamma = 1.4 # ratio of specific heats
    R = 287.052874 # gas constant for air, J/(kg K)
    rho_infty = 1.1766125010126185 # (kg/m3)
    T_infty = 300.0 # if you took very warm air then expanded it to current mach
    a = np.sqrt(gamma * R * T_infty) # speed of sound (m/s)
    # E_inlet = 1.0/(rho * (gamma-1)) + 0.5 * mach**2 
    # E_inlet = 1.0/(gamma * (gamma-1)) + 0.5 * mach**2 # cizmas in the textbook 
    cv = 1006.0 - R # J/kg k
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
    # ChatGPT says 0.25
    nu_2 = 0.0

    # 4TH ORDER DISSIPATION COEFFICIENT
    #nu4 - cizmas says 0.001 or 0.01, Jameson says 1/256, which is 2**-8
    # nu_4 = 0.00390625 # or 0.01 or 0.001
    # TAs say 0.01 and ChatGPT 4o says 0.02.
    # 0.02 works for Mach 0.3 at 11x51 grid size
    nu_4 = 0.00

    # limiting values (so that one evil cell doesn't come and kill me)
    rho_ceil = 3 # kg/m3 (reasonable)
    rho_floor = 0.00170221 # density of the air basically in space
    V_ceil = 2*a # guaranteed to not go mach 2
    V_floor = 1 # m/s - v can do whatever the firetruck it likes in lower limit
    E_ceil = 5.0e5 # a limit of about 700K (reasonable)
    E_floor = 0.5 * V_ceil**2 + 10.0*(1006.0 - R) # limit of aboutu 10K (very cold) when going max speed
    E_floor/= 3 # because I set the velocity ceiling really high and we shouldn't hit the floor anyway
    print(f'bounds set: V<{V_ceil:.1f}m/s\n\
          E>{E_floor:.2E} J/kg')

    cell_data = np.array([[[rho_infty,
                          rho_infty*(mach*np.sqrt(gamma*R*T_infty)),
                          0.0,
                          rho_infty*(cv*T_infty + 0.5*((mach*np.sqrt(gamma*R*T_infty))**2))] for j in range(len(x_list[0])-1)] for i in range(len(x_list)-1)])
    # us = np.divide(cell_data[:, :, 1].ravel(), cell_data[:, :, 0].ravel()) # get u velocities
    # vs = np.divide(cell_data[:, :, 2].ravel(), cell_data[:, :, 0].ravel()) # get v velocities
    # plot_quivers(x_list, y_list, us, vs, verbose=True)

    print(np.shape(cell_data))

    dissipation_data = np.zeros((len(x_list)-1, len(x_list[0])-1), dtype=float)
    print(f'cell data initialized with shape {np.shape(cell_data)}')
    # each index has one subtracted because for n nodes there are n-1 cells

    # state getters
    def get_rho(j,i):
        # rho
        return float64(cell_data[j,i,0])
    def get_u(j,i):
        # rho*u / rho
        return float64(cell_data[j,i,1] / cell_data[j,i,0])
    def get_v(j,i):
        # rho*v / rho
        return float64(cell_data[j,i,2] / cell_data[j,i,0])
    def get_E(j,i):
        # rho*E / rho
        return float64(cell_data[j,i,3] / cell_data[j,i,0])
    def get_T(j, i):
        """Get temperature of cell given j,i
        confirmed to agree with other calculations."""
        V = cell_data[j,i,1:3] / cell_data[j,i, 0] # vector velocity
        E = get_E(j,i)

        # get static energy
        e = E - 0.5*np.dot(V, V) # units of J/kg
        cv = 1006.0 - R # units of J/(kg K)
        T = e/cv # units of kelvin
        return np.float64(T)
    
    # nondimensional
    def get_p(j,i):
        # return ((get_rho(j,i)) * get_T(j,i)) / (T_ref) # calculated directly from other nondimensional values
        return np.float64(get_rho(j,i) * R * get_T(j,i)) # calculated by getting dimensional pressure (rho_d * R * T_d) and nondimensionalizeing as recommended by cizmas (p / rho u_ref**2)
    
    def get_H(j,i):
        # stagnation enthalpy
        # stagnation energy plus p/rho
        return np.float64(get_E(j,i) + get_p(j,i) / get_rho(j,i))

    it_number = 0
    CFL = 1.0
    while True: 
        max_res = np.zeros(4) # maximum residual in the whole iteration (numerically identify divergence of a few cells)
        ag_res = 0.0 # average residual for the current iteration
        print(it_number, end='      ')

        # dissipation calculators
        def ind(k):
            """indexify an index"""
            assert np.isclose(k%1, 0), "Index issue - irl your indeces can only be integers"
            return int(np.round(k))

        def lamb(j, i, og=-1/2):
            """ becomes lamb_xi or lamb_eta depending on what axis the off-integer value is in \n
                off integer in i  -  lamb xi  -  north velocity
                off integer in j  -  lamb eta -  east velocity
            
            input: off-integer in EITHER j or i \n
            returns: the absolute wall-normal velocity component, with respect to the wall corresponding to the half-index """


            # og determines whether to get velocity from left or right cell, negative meaning left cell, positive meaning right cell

            # bottom/left node (aka, truncated)
            j_base = int(np.round(j-epsilon)) # 5
            i_base = int(np.round(i-epsilon)) # 4

            # whether to round up or not - always round up the one with the half index
            j_round = int(np.round(j-j_base+epsilon)) # 1        0 or 1
            i_round = int(np.round(i-i_base+epsilon)) # 0           

            j_split = int(not np.isclose(i_round, 0))
            i_split = int(not np.isclose(j_round, 0))

            # node indeces
            j1 = j_base + j_round
            j2 = j_base + j_round + j_split

            i1 = i_base + i_round
            i2 = i_base + i_round + i_split

            # cell index
            #    - move up or down depending on og
            #    - move i or j depending on who is off-integer
            jc = ind(j + og*j_round )
            ic = ind(i + og*i_round)
            V = cell_data[jc, ic, 1:3] / cell_data[jc, ic, 0] # vector velocity for right/top cell

            delta_y = (y_list[j2, i2] - y_list[j1, i1])
            delta_x = (x_list[j2, i2] - x_list[j1, i1])

            # sign of the normal does not matter, only the direction
            n = np.array([delta_y, -delta_x]) / (delta_y**2 + delta_x**2)**0.5 # wall's unit normal

            a_cell =np.sqrt(gamma*R*get_T(jc, ic))

            vn_abs = np.abs(np.dot(n, V))

            if (np.isclose(j, 2) and np.isclose(i, 30.5) and np.isclose(og, 1/2)):
                print(f'taking velocity data from cell j={jc}, i={ic}\n\
                        Velocity: {V} \n\
                        Normal of {j,i}: {n}\n\
                        value of {vn_abs + a_cell}')

            # if j >= 1 and it_number%4 == 0:
            #     print(f'EIGENVALUE requested at {j,i},', end=' ')
            #     print(f'measured between the nodes (j,i = {j2, i2} and (j,i = {j1, i1})')
            #     print(f'delta x, delta_y = {delta_x, delta_y}')
            #     print(f'took velocity from cell {jc,ic}')
            #     print(f'normal: {n} with angle of {np.rad2deg(np.atan2(n[1], n[0]))}')
            #     print(f'vn_abs: {vn_abs}\n\
            #         a: {a_cell}')
            #     print(f'lambda: {vn_abs+a_cell}')
            #     plt.scatter((x_list[jc, ic]+0.05), (y_list[jc, ic]+0.05))
            #     plt.quiver(0.5*(x_list[j2, i2]+x_list[j1, i1]), 0.5*(y_list[j2, i2]+y_list[j1, i1]), n[0], n[1])
            #     plt.quiver(0.5*(x_list[j2, i2]+x_list[j1, i1]), 0.5*(y_list[j2, i2]+y_list[j1, i1]), V[0], V[1], color='r')
            #     mesh_plot(x_list, y_list, 'u')
            return vn_abs + a_cell

        def l(j,i):
            """input: off-integer indeces on either i or j \n
            returns: the length of the face using pythag"""
            j_ll = int(np.round(j-epsilon))
            i_ll = int(np.round(i-epsilon))

            delta_y = (y_list[ind(j_ll+1), ind(i_ll+1)] - y_list[ind(j_ll+int(np.round(j-j_ll+epsilon))), ind(i_ll+int(np.round(i-i_ll+epsilon)))])
            delta_x = (x_list[ind(j_ll+1), ind(i_ll+1)] - x_list[ind(j_ll+int(np.round(j-j_ll+epsilon))), ind(i_ll+int(np.round(i-i_ll+epsilon)))])

            # if it_number%4 == 0: 
            #     print(f'LENGTH (j,i={j,i}) by comparing the nodes (ji={ind(j_ll+1), ind(i_ll+1)}) and {ind(j_ll+int(np.round(j-j_ll+epsilon))), ind(i_ll+int(np.round(i-i_ll+epsilon)))} is {np.sqrt(delta_y**2 + delta_x**2)}')
            #     plt.scatter((x_list[ind(j_ll+1), ind(i_ll+1)], x_list[ind(j_ll+int(np.round(j-j_ll+epsilon))), ind(i_ll+int(np.round(i-i_ll+epsilon)))]), (y_list[ind(j_ll+1), ind(i_ll+1)], y_list[ind(j_ll+int(np.round(j-j_ll+epsilon))), ind(i_ll+int(np.round(i-i_ll+epsilon)))]))
            #     mesh_plot(x_list, y_list, np.sqrt(delta_y**2 + delta_x**2))

            return np.sqrt(delta_y**2 + delta_x**2)

        def diff_xi(obj, j, i):
            return obj(j, i+1/2) - obj(j, i-1/2)
        
        def diff_eta(obj, j, i):
            return obj(j+1/2, i) - obj(j-1/2, i)
        
        def diff2_xi(obj, j, i):
            return obj(j, i+1) - 2*obj(j, i) + obj(j, i-1)
        
        def diff2_eta(obj, j, i):
            return obj(j+1, i) - 2*obj(j, i) + obj(j-1, i)  
        
        def diff3_eta(obj, j, i):
            return diff2_eta(obj, j+1/2, i) - diff2_eta(obj, j-1/2, i)
            
        def diff3_xi(obj, j, i):
                # "at i+1/2"
                # TODO: cut out the function call for better speed
                return diff2_xi(obj, j, i+1/2) - diff2_xi(obj, j, i-1/2)

        def q(j, i):
            """input: on-integer indeces\n
                returns: cell data"""
            return cell_data[ind(j), ind(i), :]

        def switch2_eta(j,i):
            """input: off integer on j, on-integer on i \n
            in eta we vary the j"""
            # need to be evaluated "at the faces" so we need to average the left and right cell faces and average
            j_left = int(np.round(j-epsilon))
            j_right = int(np.round(j+epsilon))

            # left part of the average
            num = nu_2 * abs(diff2_eta(get_p, j_left, i))
            den = get_p(j_left, ind(i+1)) + 2*get_p(j_left, ind(i)) + get_p(j_left, ind(i-1))
            left_sw = num/den

            # right part of the average
            num = nu_2 * abs(diff2_eta(get_p, j_right, i))
            den = get_p(j_right, ind(i+1)) + 2*get_p(j_right, ind(i)) + get_p(j_right, ind(i-1))
            right_sw = num/den

            # return 0.5 * np.max([left_sw, right_sw]) # from Jameson paper
            return 0.5 * (left_sw + right_sw)

        def switch2_xi(j,i):
            """input: on integer on j, off-integer on i \n
            in xi we vary the i\n
            
            it should work fine with on-axis as well, just averages two of the same thing"""
            # need to be evaluated "at the faces" so we need to average the left and right cell faces and average
            i_left = int(np.round(i-0.1)) # i - 1/2 --> i
            i_right = int(np.round(i+0.1))# i + 1/2 --> i+1

            # override this for safety
            j_left = j
            j_right = j

            # left part of the average
            num = nu_2 * abs(diff2_xi(get_p, j_left, i_left))
            den = get_p(ind(j_left), ind(i_left+1)) + 2*get_p(ind(j_left), ind(i_left)) + get_p(ind(j_left), ind(i_left-1))
            left_sw = num/den

            # right part of the average
            num = nu_2 * abs(diff2_xi(get_p, j_right, i_right))
            den = get_p(ind(j_right), ind(i_right+1)) + 2*get_p(ind(j_right), ind(i_right)) + get_p(ind(j_right), ind(i_right-1))
            right_sw = num/den

            # return 0.5 * np.max([left_sw, right_sw]) # from jameson paper
            return 0.5 * (left_sw + right_sw)

        def switch4_xi(j,i):
            return np.max([0.0, (nu_4)-switch2_xi(j,i+1/2)])
        
        def switch4_eta(j,i):
            return np.max([0.0, (nu_4)-switch2_eta(j+1/2,i)])


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
                q_bc = 2*cell_data[j,-3, :3] - cell_data[j,-4, :3]
                cell_data[j, -2:, :3] = q_bc

                # energy
                # calculate supersonic or not
                V = cell_data[j, -2, 1:3]/cell_data[j, -2, 0]
                V_mag = np.sqrt(np.dot(V, V))

                if V_mag > a and False:
                    print(f'supersonic outlet imposed')
                    rho_E_out = get_p(j, -3) / (gamma - 1) + 0.5*get_rho(j,-2)*V_mag**2
                else:
                    rho_E_out = p_infty / (gamma - 1) + 0.5*cell_data[j, -2, 0]*V_mag**2
                
                cell_data[j,len(cell_data[0])-2:, 3] = rho_E_out
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

                        # get the momentum for the referenced cell, and reflect the normal component and preserve the tangent compoonent
                        reference_P = cell_data[domain_j, i, 1:3]
                        reflected_P = t_wall * (np.dot(reference_P, t_wall)) - n_wall * (np.dot(reference_P, n_wall))
                        cell_data[ghost_j, i, :] = np.array([cell_data[domain_j, i, 0],
                                                            reflected_P[0],
                                                                reflected_P[1],
                                                                cell_data[domain_j, i, 3]])
                        # resultant = 0.5*(cell_data[ghost_j, i, 1:3] + cell_data[domain_j, i, 1:3])
                        # print(f't_wall: {t_wall}, angle of t_wall:{np.rad2deg(np.atan2(t_wall[1], t_wall[0]))} n_wall: {n_wall}, angle of n_wall: {np.rad2deg(np.atan2(n_wall[1], n_wall[0]))}')
                        # input(f'net v at the {ghost_j}<->{domain_j} boundary at i={i}: {resultant}, angle of {np.rad2deg(np.atan2(resultant[1], resultant[0]))}')
        enforce_BCs()
  
        # display results
        # plot_scalars(x_list, y_list, cell_data[:, :, 3], caption='energy', verbose=True) # plot densities 
        # us = np.divide(cell_data[:, :, 1].ravel(), cell_data[:, :, 0].ravel()) # get u velocities
        # vs = np.divide(cell_data[:, :, 2].ravel(), cell_data[:, :, 0].ravel()) # get v velocities
        # plot_quivers(x_list, y_list, us, vs, verbose=True)
        # print out the lower boundary cells so that I can manually CFD
        if loud:
            for i in range(len(cell_data[0, :])):
                print(f'{0,i} - [', end='')
                for k in cell_data[0, i, :]:
                    print(f'{np.format_float_scientific(k, precision=3)}  ', end='')
                print(']')


        # =============== INNNER EQ =========================

        # note that this is CELL i and j not node i and j

        # len-2 because exclude right ghost, 
        # start at 2: because exclude left ghost
        # len-1 again because n-1 cells for n nodes
        for j in range(2, len(x_list)-2 - 1):
            for i in range(2, len(x_list[0])-1 - 2): 
                # print(f'\n\n\nj,i={j,i}\n\n\n')

                # calculate the area of cell (j,i) 
                A = 0.5*((x_list[j+1, i+1] - x_list[j, i])*(y_list[j+1, i] - y_list[j, i+1])  -  (y_list[j+1,i+1] - y_list[j,i])*(x_list[j+1,i] - x_list[j,i+1])) # area of the cell
                
                # calculate the fluxes between cells with averaging between the cell values

                enws = [(0, 1/2), (1/2, 0), (0, -1/2), (-1/2, 0)] # calculate f and g fluxes for east, north, west, south (this is j,i)
                for k in range(len(enws)):         
                    """ becomes lamb_xi or lamb_eta depending on what axis the off-integer value is in \n
                        off integer in i  -  lamb xi  -  north velocity
                        off integer in j  -  lamb eta -  east velocity
                    
                    input: off-integer in EITHER j or i \n
                    returns: the absolute wall-normal velocity component, with respect to the wall corresponding to the half-index """

                    j_f, i_f = j+enws[k][0], i+enws[k][1] # create boundary (off-integer) indeces
                    # print(f'k={k}, calculating flux on {j_f,i_f}')

                    # bottom/left node
                    j_base = int(np.round(j_f-epsilon)) # 5
                    i_base = int(np.round(i_f-epsilon)) # 4

                    # whether to round up or not
                    j_round = int(np.round(j_f-j_base+epsilon)) # 1        0 or 1
                    i_round = int(np.round(i_f-i_base+epsilon)) # 0           

                    # print(f'rounding bools: {j_round, i_round}')

                    # j_split = int(not np.isclose(i_round, 0))
                    # i_split = int(not np.isclose(j_round, 0))

                    # node indeces
                    # j1 = j_base + j_round
                    # j2 = j_base + j_round + j_split

                    # i1 = i_base + i_round
                    # i2 = i_base + i_round + i_split

                    # cell indeces
                    jc1 = j_base + j_round
                    ic1 = i_base + i_round

                    jc2 = j_base
                    ic2 = i_base

                    # print(f'averaging "fluxes" between {jc1, ic1} and {jc2, ic2}')

                    # calculate the f and g at the cells
                    qc = cell_data[jc1, ic1, :]
                    V = qc[1:3]/qc[0] # cell velocity
                    f1 = np.array([qc[1],
                                   qc[1]**2/qc[0] + (qc[3] - 0.5*qc[0]*np.dot(V, V))*(gamma-1),
                                   qc[1]*qc[2]/qc[0],
                                   (qc[3] + (qc[3] - 0.5*qc[0]*np.dot(V, V))*(gamma-1))*(qc[1]/qc[0])])
                    g1 = np.array([qc[2],
                                   qc[1]*qc[2]/qc[0],
                                   qc[2]**2/qc[0] + (qc[3] - 0.5*qc[0]*np.dot(V, V))*(gamma-1),
                                   (qc[3] + (qc[3] - 0.5*qc[0]*np.dot(V, V))*(gamma-1))*(qc[2]/qc[0])])

                    if loud:
                        print(f'qc({jc1,ic1}):\
                            {qc}')
                        print(f'f1: {f1}\n\
                            g1: {g1}\n')
                        

                    qc = cell_data[jc2, ic2, :]
                    V = qc[1:3]/qc[0] # cell velocity
                    f2 = np.array([qc[1],
                                   qc[1]**2/qc[0] + (qc[3] - 0.5*qc[0]*np.dot(V, V))*(gamma-1),
                                   qc[1]*qc[2]/qc[0],
                                   (qc[3] + (qc[3] - 0.5*qc[0]*np.dot(V, V))*(gamma-1))*(qc[1]/qc[0])])
                    g2 = np.array([qc[2],
                                   qc[1]*qc[2]/qc[0],
                                   qc[2]**2/qc[0] + (qc[3] - 0.5*qc[0]*np.dot(V, V))*(gamma-1),
                                   (qc[3] + (qc[3] - 0.5*qc[0]*np.dot(V, V))*(gamma-1))*(qc[2]/qc[0])])

                    if loud:
                        print(f'qc({jc2,ic2}):\
                            {qc}')
                        print(f'f1: {f2}\n\
                            g1: {g2}\n')

                    # average the vectors to get the f and g at the faces
                    # f and g are lists of the fluxes through each of the boundaries
                    f[k, :] = 0.5*(f1[:] + f2[:])
                    g[k, :] = 0.5*(g1[:] + g2[:])


                # deltas (counterclockwise)
                del_y_east = (y_list[j+1,i+1] - y_list[j,i+1]);       del_x_east = (x_list[j+1,i+1] - x_list[j,i+1])   # delta x,y in east
                del_y_north = y_list[j+1,i] - y_list[j+1,i+1];      del_x_north = (x_list[j+1,i] - x_list[j+1,i+1]) # delta x,y in north
                del_y_west = (y_list[j,i] - y_list[j+1,i]);           del_x_west = (x_list[j,i] - x_list[j+1,i]  )         # delta x,y in west
                del_y_south = y_list[j,i+1] - y_list[j,i];          del_x_south = (x_list[j,i+1] - x_list[j,i])         # delta x,y in south

                # f_res_backup = f[0]*del_y_east + f[1]*del_y_north + f[2]*del_y_west + f[3]*del_y_south <-- but faster
                f_res = f.T @ np.array([del_y_east, del_y_north, del_y_west, del_y_south])
                g_res =  g.T @ np.array([del_x_east, del_x_north, del_x_west, del_x_south])

                res = f_res - g_res

                # use JST dissipation term to reduce wiggles - https://i.imgur.com/pdUVsjX.png
                def t1(jc,ic):
                    """input - off-index values"""
                    if np.isclose(nu_2, 0.0):
                        return 0.0
                    else:
                        return switch2_xi(jc,ic) * l(jc,ic) * lamb(jc,ic,og=(i-ic) ) * diff_xi(q,jc,ic)
                def t2(jc,ic):
                    if np.isclose(nu_2, 0):
                        return 0.0
                    else:
                        return switch2_eta(jc,ic) * l(jc,ic) * lamb(jc,ic, og=(j-jc)) * diff_eta(q,jc,ic)
                def t3(jc, ic):
                    s4 = switch4_xi(jc,ic)
                    if np.isclose(s4, 0, rtol=1e-10):
                        return 0.0
                    # elif np.isclose(jc, 1.5) or np.isclose(jc, np.shape(cell_data)[0]-1.5):
                    #     return 0.0 # lambda SHOULD be zero here.
                    else:
                        return s4 * l(jc,ic) * lamb(jc, ic, og=(i-ic)) * diff3_xi(q, jc,ic)
                    
                def t4(jc, ic):
                    s4 = switch4_eta(jc,ic)
                    if np.isclose(s4, 0, rtol=1e-10):
                        return 0.0
                    # elif np.isclose(j, 1.5) or np.isclose(jc, np.shape(cell_data)[0]-1.5):
                    #     return 0.0 # lambda SHOULD be zero here.
                    else:
                        return s4 * l(jc,ic) * lamb(jc, ic, og=(j-jc)) * diff3_eta(q, jc, ic)

                try:
                    D[:] = diff_xi(t1, j,i) + diff_eta(t2, j,i) - diff_xi(t3, j, i) - diff_eta(t4, j, i)
                    dissipation_data[j,i] = D[3] 
                except IndexError:
                    print(f'found an index error at j,i={j,i}')
                    exit()

                # RK4 integration parameters
                alphas = [1/4, 1/3, 1/2, 1.0]
                # BUG: needs to be a *2 in here.
                sum_l_lamb = np.sum(np.array([lamb(jt,it, og=((j-jt)+(i-it)))*l(jt,it) for jt,it in zip([j+1/2, j, j-1/2, j],[i, i+1/2, i, i-1/2])]))

                # prevent nan residuals
                assert not np.isnan(np.dot(res, res)), "residual cannot be nan"
                assert A>0.0

                if i == 31 and j == 2: # big TE cell
                    print(f'alphas[it_number]*2*CFL / sum_l_lamb): {alphas[it_number]}*2*{CFL} / {sum_l_lamb})')
                    # something's up
                    print(f'cell lambdas: {[(name, lamb(jt,it, og=((j-jt)+(i-it)))) for name, jt,it in zip(['north ', 'east ', 'south', 'west '], [j+1/2, j, j-1/2, j], [i, i+1/2, i, i-1/2])]}')
                    for m in range(4):
                        print(f'new_q[{m}] = {cell_data[j,i,m]} - {(alphas[it_number]*2*CFL / sum_l_lamb)}*{res[m]} - {D[m]}')

                # calculate new state
                # new_q[:] = cell_data[j,i,:] - (alphas[it_number%4]*delta_t / A) * (res - D) # RK iteration with CFL
                new_q[:] = cell_data[j,i,:] - ((alphas[it_number]*2*CFL / sum_l_lamb)) * (res - D) # RK iteration without CFL
                ag_res += abs(np.linalg.norm(res))
                if max(res) >= max(max_res):
                    max_res = res

                cell_data[j,i,:] = new_q[:]
                
        # capping
        for j in range(np.shape(cell_data)[0]):
            for i in range(np.shape(cell_data)[1]):
                # cap the density
                rho = cell_data[j,i,0]
                if rho_floor > rho:# lower than the floor
                    cell_data[j,i,0] = rho_floor
                    print('WARN: rho hit floor')
                if rho_ceil < rho:  # higher than the ceiling
                    cell_data[j,i,0] = rho_ceil
                    print('WARN: capped rho')
                
                # cap the velocity and preserve direction
                V = np.sqrt(np.dot(cell_data[j,i,1:3], cell_data[j,i,1:3]))/cell_data[j,i,0]
                if V_floor > V: # less than floor
                    cell_data[j,i,1:3] = (cell_data[j,i,1:3] / V) * V_floor # normalize, then scale down to V_floor
                    print(f'WARN: velocity hit floor')
                if V_ceil < V: # higher than the ceiling
                    cell_data[j,i,1:3] = (cell_data[j,i,1:3] / V) * V_ceil # normalize, then scale up to V_ceil
                    print(f'WARN: capped velocity')
                # cap the energy
                E = cell_data[j,i,3]/cell_data[j,i,0]
                if E_floor > E:
                    cell_data[j,i,3] = E_floor*cell_data[j,i,0]
                    print(f'WARN: energy hit floor')
                if E_ceil < E:
                    print(f'WARN: capped energy')
                    cell_data[j,i,3] = E_floor*cell_data[j,i,0]

        # status updates
        # ag_res_list.append(ag_res/(len(x_list)*len(x_list[0])))
        sonic = np.max(cell_data[:, :, 1]) > a
       #  print(f'Residual Magnitude={ag_res_list[-1]:.3E}\t  max residual={np.sqrt(np.dot(max_res, max_res)):.2E}\t  sonic: {sonic}')

        # major status updates (every 8)
        if int(it_number%8) == 0 or True:

            # plt.plot(np.linspace(0, len(ag_res_list), num=len(ag_res_list)), ag_res_list)
            # plt.yscale('log');  plt.title(f'Residual Magnitude')
            # plt.savefig(f'current_residuals - jmax,imax={len(x_list)}x{len(x_list[0])}, M={mach}.png'); plt.clf()

            plot_scalars(x_list, y_list, cell_data[:, :, 0], caption=f'density - jmax,imax={len(x_list)}x{len(x_list[0])}, M={mach}') # plot densities (debugging)
            plot_scalars(x_list, y_list, cell_data[:, :, 2], caption=f'F2 - jmax,imax={len(x_list)}x{len(x_list[0])}, M={mach}') # plot densities (debugging)

            plot_scalars(x_list, y_list, np.divide(cell_data[:, :, 3], cell_data[:, :, 0]), caption=f'energy - jmax,imax={len(x_list)}x{len(x_list[0])}, M={mach}') # plot densities (debugging)
            plot_scalars(x_list, y_list, dissipation_data[:, :], caption=f'dissipations - jmax,imax={len(x_list)}x{len(x_list[0])}, M={mach}') # plot densities (debugging)

            # x velocities and y velocities
            us = np.divide(cell_data[:, :, 1].ravel(), cell_data[:, :, 0].ravel()) # get u velocities
            vs = np.divide(cell_data[:, :, 2].ravel(), cell_data[:, :, 0].ravel()) # get v velocities
            plot_quivers(x_list, y_list, us, vs, caption=f'velocities - jmax,imax={len(x_list)}x{len(x_list[0])}, M={mach}', verbose=True)
            machs = np.power(np.power(us, 2) + np.power(vs, 2), 0.5) / a
            plot_scalars(x_list, y_list, machs, caption=f'Machs - jmax,imax={len(x_list)}x{len(x_list[0])}, M={mach}')

            # ravelled_size = np.shape(cell_data)[0] * np.shape(cell_data)[1] * np.shape(cell_data)[2]
            # longer_size = int(max(ravelled_size, len(ag_res_list)))
            
            # data_padded = cell_data.ravel().tolist() + [np.nan] * (longer_size - ravelled_size)
            # residuals_padded = ag_res_list + [np.nan] * (longer_size - len(ag_res_list))

            # df = pd.DataFrame({'data': data_padded,
            #                    'convergence': residuals_padded})
            # df.to_csv(f'j_max={len(x_list)-1},i_max={len(x_list[0])-1},M={round(mach, 1)}.csv')
            # del df # aha!
            # input()
        it_number += 1 # my iteration number just for my own bookkeeping

    exit() # only 1 iteration

    # print out the force on the bottom bump as well as a mach plot
    corner_pts = [find_arg_nearest(x_list[0], 2.0), find_arg_nearest(x_list[0], 3.0)]
    Fx = 0.0
    Fy = 0.0
    for i in range(corner_pts[0], corner_pts[1]):
        j = 2
        delta_y = (y_list[ind(j+1), ind(i+1)] - y_list[ind(j+int(np.round(j-j+epsilon))), ind(i+int(np.round(i-i+epsilon)))])
        delta_x = (x_list[ind(j+1), ind(i+1)] - x_list[ind(j+int(np.round(j-j+epsilon))), ind(i+int(np.round(i-i+epsilon)))])
        n = np.array([-delta_y, delta_x]) / (delta_y**2 + delta_x**2)**0.5
        F = get_p(j,i) * l(j-0.5,i) * n
        Fx += F[0]
        Fy += F[1]
    msg = f'X force: {Fx} N\n\
            Y force: {Fy} N'
    print(msg)
    with open(f'magic number- M={mach}, ji={len(x_list),len(x_list[0])}.txt', 'w') as my_glorious_result:
        my_glorious_result.write(msg)

    return 0


# main sequence
def main():

    machs = [0.3]
    for M in machs:
        print(f'== MACH {M} ==')
        shapes = [(10,50)]

        for shape in shapes:
            # get the mesh
            mesh_fname = f'mesh_sh={shape[0]}x{shape[1]}.csv'

            # generate mesh
            my_mesh = Mesh(shape)
            try:
                my_mesh.read(mesh_fname)
            except FileNotFoundError:
                my_mesh.generate()
                my_mesh.ghostify() # saved mesh must be ghostified.
                my_mesh.save(mesh_fname)
                my_mesh.plot() 
            
            df = pd.read_csv(mesh_fname)
            x_list = my_mesh.x_list
            y_list = my_mesh.y_list
            
            # run the simulation and generate results
            run(x_list, y_list, M)

            # TODO: hand calc if the density goes up or down at the choke point for mach 0.3

    print(f'Process finished with exit code 0')
    exit()

# run the script
if __name__ == "__main__":
    main()