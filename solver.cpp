#include <iostream> // console messages
#include <vector> // vectors my beloved
#include "csv.hpp"
#include "solver.h"
#include <cmath>
// TODO: rename to solve_engine or something like that

// data inputter
void Solution::innit(arrayD3 _mesh_data) {

    /*there is a bug in here that is making the code put the wrong amount of states per cell*/

    std::cout << "innit\n";

    // take the mesh data in
    mesh_data = _mesh_data;

    std::cout << "j_max, i_max = " << j_max << "," << i_max << "\n";
    system("pause");

    q.resize(j_max);
    // fill in the initial state q (go up to the max node index MINUS ONE because cells are one less than nodes)
    for (int j = 0; j<=j_max-1; j++) { 
        q[j].resize(i_max);
        for (int i = 0; i<=i_max-1; i++) {
                // DEBUGGING: print out indeces to ensure iteration is correct
                std::cout << "current index (j,i) = " << j << "," << i << "\n";
                q[j][i] = std::vector<float>{
                    static_cast<float>(p_infty / (R * t_infty)), // solve for density
                    static_cast<float>((p_infty / (R * t_infty)) * (mach_infty * sqrt(1.4*R*t_infty))),
                    0.0f,
                    static_cast<float>(0.5*pow(mach_infty * sqrt(1.4*R*t_infty), 2) + (cv*t_infty))
                };
        }
    }

    // DEBUGGING: print out the mesh_data to see if it is good
    for (const auto& row : q) {
        for (const auto& state: row) {
            // for (const auto& xy: pair) {
            //     std::cout << xy;
            // }
            std::cout << state[0] << " ";
        }
        std::cout << "\n";
    }
}

// simple state getters
float Solution::p(int j, int i) {
    std::cout << "Pressure getter reporting";
    return 0.0;
}                 
float Solution::T(int j, int i) {
    std::cout << "temperature getter reporting";
    return 0.0;
}                 
float Solution::rho(int j, int i) {
    std::cout << "density getter reporting";
    return 0.0;
}
float Solution::e(int j, int i) {
    std::cout << "specific internal static energy getter reporting";
    return 0.0;
}


// internal functions
float Solution::l(float j, float i) {
    // returns the length of a cell wall given the wall's index in off-integer notation
    std::cout << "length function reporting";
    // convert to wall-coordinates
    return 0.0;
    // returns the length of a cell wall  
}

float Solution::lambda(float j, float i) {
    /* input: j and i in off-integer notation
     body:
         takes the velocity AT THE WALL (average) between the two cells
         takes l = the wall normal of the velocity
     returns:
         l
    */
    std::cout << "lambda function reporting";
    // returns the eigenvalue at the cell
    return 0.0;
}

float Solution::D(int j, int i) {
    float nu_2 {0.25};
    float nu_4 {0.002};
    /*input: 
        j and i cells
    body;
        calculates second order dissipations
        calculates fourth order dissipations
    returns:
        total dissipations*/
    
    std::cout << "dissipation function reporting";
    return 0.0;
}

void Solution::boundary_conditions() {
    /*input:
        None
    body: 
        edits q such that boundary conditions are enforced
    output: 
        None*/
    std::cout << "boundary condition enforcement function reporting\n";
}

arrayD3 Solution::get_q() {
    /*inputs:
        none
    outputs:
        q   */
    return q;
}

void Solution::iterate() {
    /* conduct one iteration*/
    std::cout << "Iterate\n";

    /*
    boundary_conditions();
    for (int j = 0; j<j_max; j++) {
        for (int i = 0; i<i_max; i++) {

        }
    }
    */

}
