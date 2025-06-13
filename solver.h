#ifndef SOLVER_H
#define SOLVER_H
#include <vector>

// macros

// float array of dimension 3 columns of rows of (x,y)
// float array of dimension 3: columns of rows of (rho, rho*u, rho*v, rho*E)
using arrayD3 = std::vector<std::vector<std::vector<float>>>;
using arrayD2 = std::vector<std::vector<float>>;

#define CFL 1.0
#define R 287.052
#define pi 3.1415
#define gamma 1.4
#define cp (gamma*R/(gamma-1))
#define cv (cp-R)
// header - definitions go here for clarity

// setup.cpp
void get_config();
arrayD3 get_mesh_data();
void save_data(arrayD3);

// can't declare variables here, or else it appears in all the files, leading to a conflict.
// extern is supposed to fix this definition/declaration confusion?
extern int i_max; // max index of the ghost nodes
extern int j_max; // max index of the ghost nodes
extern int max_iterations; // maximum number of iterations to perform
extern float t_infty; // kelvin
extern float p_infty; // pascals
extern float mach_infty;
extern std::string mesh_fname;
extern std::string res_fname;

// solver.cpp
class Solution {
    public:
        // variables
        float ag_res;                                       // aggregate residual
        int iteration_count;

        void innit(arrayD3 _mesh_data);            // take in mesh data and initialize the state
        void iterate();                                     // do 1 time-step
        arrayD3 get_q();                                    // return the results

    private:
        // private data
        arrayD3 mesh_data;                      // private - do not touch my data.
        arrayD3 q;                              // 3D state array
        arrayD3 new_q;                 // 3D state array for intermediate (RK45 iterations)
        arrayD3 f;
        arrayD3 g;
        
        // private functions (internal workings)

        void update_BCs();               // enforce boundary conditions
        void update_f_g();

        // simple state getters
        float p(int i, int j);                  // get pressure at a cell
        float T(int i, int j);                  // get temperature at a cell
        float rho(int i, int j);                // get density at a cell
        float e(int i, int j);                  // get specific internal static energy at a cell
    
        // advanced numerical finders
        float l(int i, int j, float off_i, float off_j);              // length of cell wall (input off-integer values)
        float lambda(float i, float j);         // eigenvalue of cell at the wall (input off-integer values)
        float cizmas_lambda(int i, int j, float off_i, float off_j);

        std::vector<float> D(int i, int j);              // dissipation at the cell
    };

#endif