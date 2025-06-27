# ABOUT
This project started as a final project for the class AERO415. I initially wrote a CFD solver in Python, because that was what I was familiar with. Even after pouring about a week of 14+ hour, highly focused workdays into the project, it still had two issues:
 1) slowness - The Python script took multiple hours to run a satisfactory number of iterations, even in the most crude case. This was after spending a significant amount of time trying to optimize the Python. I investigated Cython, CuPy, and NumPy vectorization. I recall the script taking about twenty seconds to run a SINGLE iteration on 50x200 mesh. This is unacceptible. After recreating my code in C++, I found the unoptimized C++ to be about 300x faster than the unoptimized Python.
 2) numerical divergence - Even after the multiple hours of running, the solver was unable to converge on a result, and the residuals would eventually increase to infinity.

While a lot of people would have walked away and never looked back (since the semester was over, after all, and any grade incintive would have been eliminated), I could not walk away bested by this challenge. I chose to finish this code for the same reason that I chose to study aerospace engineering in the first place. In the words of President Kennedy, I do these things ***"not because they are easy, but because they are hard."***

# OVERALL STRUCTURE
## How to Run
1) compile the solver in build using the VSCode Cmake and CMakeLists.txt
2) run main.py, which will:
    - generate a mesh 
    - edit configuration file solver_commands.csv
    - call solver_engine.exe
    - display results


## Code Structure
Python: Main.py generates mesh and sets up config file

C++:
 - declared functions and classes in the .h file
 - in subordinate files, filled in definitions
 - subordinate files
    - main.cpp - overall logic/flow for the program
    - setup.cpp - get data from config file (solver_commands.csv), and open the mesh file (name is given in the config file)
    - solver.cpp - Houses class Solution
        - Innit() - initializes the state vector, temporary state vector, f and g fluxes, and dissipation with initial values
        - p() - calcualtes the pressure at a given cell
        - T() - calculates the temperature at a given cell
        - l() - calculates the length of a face
        - lambda() - calcualtes the cell eigenvalue (fastest speed that a wave can travel) at a cell's face
        - switches - Detects abrupt pressure changes (shocks) and changes to second order dissipation
        - update_D() - Updates dissipation at one cell 
        - Update_BCs() - Updates all boundary conditions
        - update_f_g - updates all the flux values
        - get_q() - returns a read-only value of the system's state
        - iterate() - conduct one full iteration of the solution


# ABOUT CELL INDEXING
 - i be associated with the axis 0 (which will be most closely associated with y, or the vertical direction)
 - j be associated with axis 1 (which will be most closely associated with x, or the horizontal direction)
 - k will be associated with the state index, with the indeces corresponding to the values like so:
    - 0 - density
    - 1 - volume-specific x momentum
    - 2 - volume-specific y momentum
    - 3 - volume-specific energy
## indexing example: 10x50 grid:
### normal mesh:
cells
 - 10x50 cells
 - cell indeces span 0-9, 0-49

nodes
 - 11x51 nodes
 - node indeces span 0-50

### mesh when ghost cells are added:
everything gets a +4
 - 14x54 cells
 - cell indeces span 0-13, 0-53


## converting indeces:


|     | Cell number           | Cell index               | Node number           | Node Index         |
| --- | --------------------- | ------------------------ | --------------------- | ------------------ |
| +5  |                       |                          | number of ghost nodes |                    |
| +4  | number of ghost cells |                          |                       | max index of ghost nodes |
| +3  |                       | max index of ghost cells |                       |                    |
| +2  |                       |                          |                       |                    |
| +1  |                       |                          | number of core nodes  |                    |
| +0  | number of core cells  |                          |                       | max index of nodes |
| -1  |                       | max index of core cells  |                       |                    |


## off-integer notation for cell walls
Cells:
 - Begin at 0,0 in the lower left corner
 
Cell walls:
 - Begin at -1/2, 0 for the bottom-most border in the lower left cell
 - Begin at 0, -1/2 for the left-most border in the lower left cell.
 

 # SPECIAL THANKS 

## Libraries
 - CSV reader: https://github.com/vincentlaucsb/csv-parser

 ## TO THE DEBUGGING DUCK
 - https://i.imgur.com/j3b668J.png, for listening to me explain my program


## To the TAMU CFD Lab
 - Dr. Paul Cizmas, https://engineering.tamu.edu/aerospace/profiles/pcizmas.html, for his notes
 - To Matthew Schultz, Ciprian Comsa, and Jay Standridge, for their encouragement and advice


