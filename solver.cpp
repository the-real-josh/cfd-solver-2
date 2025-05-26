#include <iostream> // console messages
#include <vector> // vectors my beloved
#include "external/csv-parser/include/csv.hpp" // csv reading (how to import??)
#include "solver.h"

// data inputter
void Solution::innit(arrayD3 _mesh_data, float M) {
    // take the mesh data in
    mesh_data = _mesh_data;

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
    // chat I think it should be the vector magnitude of the velocity not the scalar absolute value of the u

    std::cout << "lambda function reporting";
    // returns the eigenvalue at the cell
    return 0.0;
}

float Solution::D(int j, int i) {
    std::cout << "dissipation function reporting";
    return 0.0;
}

void Solution::boundary_conditions() {
    std::cout << "boundary condition enforcement function reporting";
}

void Solution::iterate() {
    boundary_conditions();

}


arrayD3 Solution::get_q() {
    // data getter
    return q;
}
