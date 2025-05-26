#ifndef SOLVER_H
#define SOLVER_H

#include <vector>

// float array of dimension 3 columns of rows of (x,y)
// float array of dimension 3: columns of rows of (rho, rho*u, rho*v, rho*E)

using arrayD3 = std::vector<std::vector<std::vector<float>>>;

// header - definitions go here for clarity

// setup.cpp
arrayD3 get_mesh_data();

// main.cpp
class Solution {
    public:
        // variables
        float ag_res;                                       // aggregate residual

        void innit(arrayD3 _mesh_data, float M);                     // take in mesh data and initialize the state
        void iterate();                                     // do 1 time-step
        arrayD3 get_q();                                           // return the results

    private:
        // private data
        arrayD3 mesh_data;                      // private - do not touch my data.
        arrayD3 q;                              // 3D state array (for debugging)
        arrayD3 new_q;
        // private functions (internal workings)

        void boundary_conditions();               // enforce boundary conditions

        // simple state getters
        float p(int j, int i);                  // get pressure at a cell
        float T(int j, int i);                  // get temperature at a cell
        float rho(int j, int i);                // get density at a cell
        float e(int j, int i);                  // get specific internal static energy at a cell
    
        // advanced numerical finders
        float l(float j, float i);              // length of cell wall (input off-integer values)
        float lambda(float j, float i);         // eigenvalue of cell at the wall (input off-integer values)
    
        float D(int j, int i);              // dissipation at the cell
    };

#endif