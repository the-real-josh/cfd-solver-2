#ifndef WHAT GOES HERE
#define WHAT GOES HERE

#include <vector>

// float array of dimension 3 columns of rows of (x,y)
// float array of dimension 3: columns of rows of (rho, rho*u, rho*v, rho*E)
#define arrayD3 std::vector<std::vector<std::vector<float>>> 

// header - definitions go here for clarity

// setup.cpp
arrayD3 get_mesh_data();
void boundary_conditions(float mach);
void initialize(float mach);

// main.cpp
class Solution {
    public:
        float ag_res;   // aggregate residual
        void iterate(); // do 1 time-step
        void save();    // save the results
        void innit();
    private:
        arrayD3 mesh_data;
        arrayD3 q;                              // 3D state array (for debugging)

        // simple state getters
        float p(int j, int i);                  // get pressure at a cell
        float T(int j, int i);                  // get temperature at a cell
        float rho(int j, int i);                // get density at a cell
        float e(int j, int i);                  // get specific internal static energy at a cell
    
        // advanced numerical finders
        float l(float j, float i);              // length of cell wall (input off-integer values)
        float lambda(float j, float i);         // eigenvalue of cell at the wall (input off-integer values)
    };
