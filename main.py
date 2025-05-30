import subprocess
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

class Mesh:
    def __init__(self, dimensions):
        self.dimensions = dimensions
    def generate(self):
        min_residual=1e-6

        """# Resolution:
            - resolution of under 4 causes everything to completely break
            - resolution is nodes per unit
            - change the resolution to prove mesh-independence of solution

            # min residual:
            - the average amount that is changed by the iterator
            - default is 1e-6
            """
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
        J_max = self.dimensions[0] # number of nodes in the y direction
        I_max = self.dimensions[1] # for debuggin'

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
        max_iterations = 100*self.dimensions[0] # iteration ceiling
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
        delta_x = np.average(np.diff(self.x_list[0, :]))
        delta_y = np.average(np.diff(self.y_list[:, 0].ravel()))

        # append in the top and bottom
        self.x_list = np.concatenate(([self.x_list[0, :]], [self.x_list[0, :]], self.x_list[:, :], [self.x_list[-1, :]], [self.x_list[-1, :]]), axis=0)
        self.y_list = np.concatenate(([self.y_list[0, :] - 2*delta_y], [self.y_list[0, :] - delta_y], self.y_list[:, :], [self.y_list[-1, :] + delta_y], [self.y_list[-1, :] + 2*delta_y]), axis=0)

        # append in the left and the right
        self.x_list = np.concatenate((self.x_list[:, 0:1] - 2*delta_x, self.x_list[:, 0:1] - delta_x, self.x_list[:], self.x_list[:,  -1:] + delta_x, self.x_list[:,  -1:] + 2*delta_x), axis=1)
        self.y_list = np.concatenate((self.y_list[:, 0:1], self.y_list[:, 0:1], self.y_list[:], self.y_list[:,  -1:], self.y_list[:,  -1:]), axis=1)

    def save(self):
        df = pd.DataFrame({f'{self.dimensions}-x':self.x_list.ravel(),
                            f'{self.dimensions}-y':self.y_list.ravel()})
        df.to_csv()

    def read(self):
        df = pd.read_csv('saved_meshes.csv')
        self.x_list = df[f'{self.dimensions}-x']
        self.y_list = df[f'{self.dimensions}-y']
        pass
    def plot(self):
        """Pass data in as x_list and y_list as a list of lists in the shape of the mesh"""
        # fancy plot with skewness
        squishedness = 0
        n = 0
        i_ubound = len(self.x_list)
        j_ubound = len(self.x_list[0])
        for i in range(i_ubound):
            for j in range(j_ubound):
                try:
                    if (i==i_ubound and j==j_ubound) or (i==i_ubound-1 and j==j_ubound-1):
                        pass
                    elif j==j_ubound-1:
                        plt.plot((self.x_list[i, j], self.x_list[i+1, j]), (self.y_list[i, j], self.y_list[i+1, j]), color='black')
                    elif i==i_ubound-1:
                        plt.plot((self.x_list[i, j], self.x_list[i, j+1]), (self.y_list[i, j], self.y_list[i, j+1]), color='black')
                    else:
                        plt.plot((self.x_list[i, j], self.x_list[i+1, j]), (self.y_list[i, j], self.y_list[i+1, j]), color='black')
                        plt.plot((self.x_list[i, j], self.x_list[i, j+1]), (self.y_list[i, j], self.y_list[i, j+1]), color='black')
                        horiz = np.array([self.x_list[i+1, j], self.y_list[i+1, j]]) - np.array([self.x_list[i, j], self.y_list[i, j]])
                        vert = np.array([self.x_list[i, j+1], self.y_list[i, j+1]]) - np.array([self.x_list[i, j], self.y_list[i, j]])
                        if not 0:
                            squishedness += np.abs(np.dot(horiz, vert) / (np.linalg.norm(horiz)*np.linalg.norm(vert)))
                            n += 1
                except IndexError:
                    print(f'grapher - index error @ji={j,i} for shape of {np.shape(self.x_list)}')
                    pass
        plt.title(f'mesh: {self.dimensions}')
        plt.scatter((2, 2, 3, 3), (0, 1, 0, 1))
        plt.axis('equal')
        plt.show()
        self.squishedness = squishedness/(i_ubound*j_ubound)


def run(M):
    # update commands.txt
    with open('commands.txt', 'w') as f:
        f.write(f'shape: 3,3\n\
                mach: {M}\n\
                iterations: {1}')
        
    # maybe later include some stuff about conditions in here
    
    # call up mr c++
    subprocess.call(["build/run.exe"])

def results():
    # plot results

    # plot density:
    

    return 0

def main():
    shapes = [(11, 21)]
    machs = [0.5]
    for shape in shapes:
        my_mesh = Mesh(shape)
        try:
            my_mesh.read()
        except:
            my_mesh.generate()
            my_mesh.ghostify() 
            my_mesh.save()
        my_mesh.plot()
        input(f'mesh acquisition complete')
        for mach in machs:
            run(mach)
        results()

if __name__ == "__main__":
    main()