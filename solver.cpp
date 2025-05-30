#include <iostream> // console messages
#include <vector> // vectors my beloved
#include "csv.hpp"
#include "solver.h"

// TODO: rename to solve_engine or something like that

// data inputter
void Solution::innit(arrayD3 _mesh_data, float M) {
    // take the mesh data in
    mesh_data = _mesh_data;
    I_max = mesh_data[0].size();
    J_max = mesh_data.size();

    // fill in the initial state q

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
    boundary_conditions();
    for (int j = 0; j<J_max; j++) {
        for (int i = 0; i<I_max; i++) {

        }
    }

}
