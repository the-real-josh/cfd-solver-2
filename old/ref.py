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


# main sequence
def main():
    # ensure directory exists
    if not os.path.exists(f'{savedir}'):
        os.makedirs(f'{savedir}')


    machs = [0.3]#, 0.3, 0.7]
    for M in machs:
        print(f'== MACH {M} ==')
        # shapes = [(11,51)]#, (21,101), (41,201)]
        shapes = [(11,51)]#, (21,101), (11, 51)]
        # shapes = [(9, 25)]

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
            
            # run the simulation and generate results
            # run(x_list, y_list, M)
            import src.solver
            src.solver.main(shape, M)

    print(f'Process finished with exit code 0')
    exit()

# run the script
if __name__ == "__main__":
    main()