# how to run:
1) compile the solver in build using the VSCode Cmake and CMakeLists.txt
2) run main.py, which will:
    - generate a mesh 
    - edit configuration file solver_commands.csv
    - call solver_engine.exe
    - display results

# about cell indexing
The cell indexing scheme in the Python mesher is made of a vector of row-containing vectors. The value at indeces (i,j) may be accessed by cell_data[j][i]

The cell indexing scheme in the C++ is made of a vector of column-containing vectors, such that the value at indeces (i,j) may be accessed q[i][j]

this will also be followed by the file naming convention

**In both**, i is associated with x, and j is associated with y.

The preferred notation will henceforth be a vector of column-containing vectors, because that is more intuitive.  

## Example: 10x50 grid:
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


## converting:


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
 