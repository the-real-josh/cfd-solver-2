import subprocess
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

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

        def find_arg_nearest(array,value):
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
        J_max = self.core_dimensions[0]+1 # y direction
        I_max = self.core_dimensions[1]+1 # x direction

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

        # preliminary calcs
        delta_xi = 1/(I_max);   delta_eta = 1/(J_max) # this is because xi and eta are defined to be on the interval {0<x<1}. In python if your list has max index of 3, then you have 2 divisions

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
        ghostified_dimensions = [i+4+1 for i in self.core_dimensions]
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


def run(mesh_core_dimensions=None, mach=None, mesh_fname=None, results_fname=None):
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

def plot_results(mesh_core_dimensions=None, mesh_fname=None, results_fname=None, verbose=True):
    """
    Plots solver results
        plots density and energy with a countourf and colorbar scale
        docs: https://www.geeksforgeeks.org/matplotlib-pyplot-contourf-in-python/

        plots velocities with a quiver
        docs: https://matplotlib.org/stable/gallery/images_contours_and_fields/quiver_simple_demo.html#sphx-glr-gallery-images-contours-and-fields-quiver-simple-demo-py
    """

    # get the pre-ghostified mesh
    curr_mesh = Mesh(mesh_core_dimensions)
    curr_mesh.read(mesh_fname)
    x = curr_mesh.x_list
    y = curr_mesh.y_list

    # calculate centers
    data_dimensions = [k+4 for k in mesh_core_dimensions] # cells have maximum indeces one less
    x_cell_centers = np.zeros(data_dimensions)
    y_cell_centers = np.zeros_like(x_cell_centers)
    for j in range(data_dimensions[0]):
        for i in range(data_dimensions[1]):
            x_cell_centers[j,i] = np.mean(np.array([x[j,i],
                                                    x[j+1,i],
                                                    x[j,i+1],
                                                    x[j+1,i+1]]))
            y_cell_centers[j,i] = np.mean(np.array([y[j,i],
                                                y[j+1,i],
                                                y[j,i+1],
                                                y[j+1,i+1]]))
            
    # test plot with no cfd data
    plt.scatter(x_cell_centers, y_cell_centers)
    mesh_plot(x,y)
    plt.show()

    # extract results
    df = pd.read_csv(results_fname)
    rho = df['rho'].to_numpy().reshape(data_dimensions)
    u = np.divide(df['rho_u'].to_numpy().reshape(data_dimensions), rho)
    v = np.divide(df['rho_v'].to_numpy().reshape(data_dimensions), rho)
    E = np.divide(df['rho_E'].to_numpy().reshape(data_dimensions), rho)

    # plot density
    mesh_plot(x,y)
    z = np.ma.masked_where(rho <= 0, rho)
    cs = plt.contourf(x_cell_centers, y_cell_centers, z,
                      locator=matplotlib.ticker.LogLocator(),
                      cmap="bone")
    cbar = plt.colorbar(cs)
    plt.title("density contour plot attempt")

    # plot density with a countourf and colorbar scale
    # https://www.geeksforgeeks.org/matplotlib-pyplot-contourf-in-python/

    # plot velocities with a quiver
    # https://matplotlib.org/stable/gallery/images_contours_and_fields/quiver_simple_demo.html#sphx-glr-gallery-images-contours-and-fields-quiver-simple-demo-py
    if verbose:
        plt.show()
    else:
        plt.savefig(mesh_fname + "density")
        plt.clf()

    # plot energy same as density

    # end
    del df
    plt.close()


# for n cells:
#  - max cell index: n-1
#  - max node index: n
#  - number of nodes: n+1
def main():
    # generate mesh
    core_dimensions = [(4, 20), (10,50)] # core mesh shapes
    machs = [0.5]
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

        my_mesh.plot()

        # run solver
        for mach in machs:
            results_fname = f'results_sh={c_dim[0]}x{c_dim[1]} M={mach}.csv'

            run(mesh_core_dimensions=c_dim, mach=mach, results_fname=results_fname, mesh_fname=mesh_fname)
        
            # view results
            plot_results(mesh_fname, results_fname)

def test_main():
    df = pd.read_csv('test.csv')
    x = df['x'].to_numpy().reshape((7,6))
    y = df['y'].to_numpy().reshape((7,6))
    
    # each column of points is a vector, and the array is vector of vectors.
    # still correct that x -> i, and y -> j

    mesh_plot(x,y); plt.show()

    run(mesh_core_dimensions=(2,1),
        mach=0.5,
        results_fname='output.csv',
        mesh_fname='test.csv')
    
    plot_results(mesh_core_dimensions=(2,1), mesh_fname='test.csv', results_fname='output.csv')

if __name__ == "__main__":
    test_main()