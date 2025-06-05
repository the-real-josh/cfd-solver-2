import subprocess
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import pyvista as pv

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


def run(mesh_core_dimensions=None, mach=None, mesh_fname=None, results_fname=None, max_iterations=100):
    """ - Edit a configuration file ("./solver_commands.csv")
        - call the solver: ./out/build/default/solver_engine.exe"""

    # other givens required to fix the state
    T_infty = 300.0
    p_infty = 101325.0

    # update commands.txt

    # NOTE about determining i_max and j_max from core dimensions.
    # +4 (see readme.txt#about cell indexing##converting)
    d = {'i_max': mesh_core_dimensions[0] + 4, 
         'j_max': mesh_core_dimensions[1] + 4,
         'max_iterations': max_iterations, 
         't_infty': T_infty,
         'p_infty': p_infty,
         'mach_infty': mach,
         'mesh_fname': mesh_fname,
         'res_fname': results_fname}
    df = pd.DataFrame([d])
    print(df)
    df.to_csv("solver_commands.csv")
    del df

    # call up mr c++
    subprocess.call(["out/build/default/solver_engine.exe"])

    return 0

def plot_results_pv(mesh_core_dimensions=None, mesh_fname=None, results_fname=None, verbose=True):
    savedir = '/'

    # -------- get data from the files -------------------------------
    # get the pre-ghostified mesh
    curr_mesh = Mesh(mesh_core_dimensions)
    curr_mesh.read(mesh_fname)
    x = curr_mesh.x_list
    y = curr_mesh.y_list

    # extract results
    df = pd.read_csv(results_fname)
    data_dimensions = [k+4 for k in mesh_core_dimensions] # cells have maximum indeces one less
    rho = df['rho'].to_numpy().reshape(data_dimensions)
    u = np.divide(df['rho_u'].to_numpy().reshape(data_dimensions), rho)
    v = np.divide(df['rho_v'].to_numpy().reshape(data_dimensions), rho)
    E = np.divide(df['rho_E'].to_numpy().reshape(data_dimensions), rho)




    # -------- quiver-plot velocities -------------------------------
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
    u = u.ravel()
    v = v.ravel()

    speeds = np.pow(np.pow(u, 2) + np.pow(v, 2), 0.5).reshape(np.shape(x)[0]-1, np.shape(y)[1]-1) # plot
    grid.cell_data['velocities'] = speeds.ravel()
    
    # Visualize the mesh
    plotter = pv.Plotter(off_screen=(not verbose))
    plotter.add_mesh(grid, scalars='velocities', cmap='plasma', show_edges=True) # cmaps = plasma, viridis
    
    centers = np.array(zip2(x_centers, y_centers, np.zeros_like(x_centers)))
    velocities = np.array(zip2(u, v, np.zeros_like(u)))
    plotter.add_arrows(centers, velocities, mag=0.1/np.max(speeds.ravel()))
    plotter.camera_position = 'xy'
    plotter.show_axes()

    # give output
    if not verbose:
        plotter.screenshot(f'velocities.png')
    elif verbose:
        plotter.show()

    plotter.close()
    del plotter




    # ------------------ Plot scalars  ------------------
    # TODO: add other scalars: temperature, pressure, machs
    captions = ['Density', 'Energy']
    datas = [rho, E]
    for data, caption in zip(datas, captions):

        # Create a 2D structured grid object
        grid = pv.StructuredGrid()
        grid.dimensions = x.shape[1], x.shape[0], 1     # Set the dimensions of the grid
        pv.set_plot_theme("document")
        # Create the points array with z coordinates as zeros
        z = np.zeros_like(x.flatten())
        points = np.stack((x.flatten(), y.flatten(), z), axis=-1)
        # Set the points of the grid
        grid.points = points
        grid.cell_data[caption] = data.ravel()

        # only put edges when the cells are reasonably large
        if len(x) > 30:
            edges = False
        else:
            edges = True

        plotter = pv.Plotter(off_screen=(not verbose))
        plotter.add_mesh(grid, scalars=caption, cmap='plasma', show_edges=edges) # default colormap is viridis
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

# for n cells:
#  - max cell index: n-1
#  - max node index: n
#  - number of nodes: n+1
def main():
    # generate mesh
    # i is no longer 
    core_dimensions = [(10,50)] # core mesh shapes
    machs = [0.3]
    for c_dim in core_dimensions:

        # mesh file name
        mesh_fname = f'mesh_sh={c_dim[0]}x{c_dim[1]}.csv'

        # generate mesh
        my_mesh = Mesh(c_dim)
        try:
            my_mesh.read(mesh_fname)
        except FileNotFoundError:
            my_mesh.generate()
            my_mesh.ghostify() # saved mesh must be ghostified.
            my_mesh.save(mesh_fname)
            my_mesh.plot()

        # my_mesh.plot() # 

        # run solver
        for mach in machs:
            results_fname = f'results_sh={c_dim[0]}x{c_dim[1]} M={mach}.csv'

            # put this back later
            run(mesh_core_dimensions=c_dim,
                 mach=mach,
                   results_fname=results_fname,
                     mesh_fname=mesh_fname,
                     max_iterations=250)
        
            # for now, only view result.
            plot_results_pv(mesh_core_dimensions=c_dim, mesh_fname=mesh_fname, results_fname=results_fname)

def test_main():
    """This test case shows some strange triangulation done by PyVista. I don't like PyVista, but it's what we have."""
    df = pd.read_csv('test.csv')
    x = df['x'].to_numpy().reshape((7,6)) # i_max = 7, j_max = 6
    y = df['y'].to_numpy().reshape((7,6))

    mesh_plot(x,y); plt.show()

    run(mesh_core_dimensions=(2,1), # core i_max = 1, core j_max = 2
        mach=0.5,
        results_fname='output.csv',
        mesh_fname='test.csv')
    
    plot_results_pv(mesh_core_dimensions=(2,1), mesh_fname='test.csv', results_fname='output.csv')

if __name__ == "__main__":
    main()