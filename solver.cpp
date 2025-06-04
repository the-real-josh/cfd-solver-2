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

    std::cout << "i_max, j_max = " << i_max << "," << j_max << "\n";
    // system("pause");

    q.resize(i_max);
    // fill in the initial state q (go up to the max node index MINUS ONE because cells are one less than nodes)
    for (int i = 0; i<=i_max-1; i++) { 
        q[i].resize(j_max);
        for (int j = 0; j<=j_max-1; j++) {
                q[i][j] = std::vector<float>{
                    static_cast<float>(p_infty / (R * t_infty)), // solve for density
                    static_cast<float>((p_infty / (R * t_infty)) * (mach_infty * sqrt(gamma*R*t_infty))),
                    0.0f,
                    static_cast<float>((p_infty/(R*t_infty))*(0.5*pow(mach_infty * sqrt(gamma*R*t_infty), 2) + (cv*t_infty)))
                };
        }
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
    std::cout << "Enforcing Boundary Conditions\n";

    // temporary variables for enforcement
    std::vector<float> bdry_velocity;
    std::vector<float> wall_vec;
    std::vector<float> wall_norm;
    float v_dot_n;
    int i_bdry; // i value of the cell that is mirrored by the ghost cell (boundary cell)
    int i_ghost;    // i value of the ghost cell

    // inlet boundary condition (compiles and runs good)
    std::vector<float> q_inlet = {
        static_cast<float>(p_infty/(R*t_infty)),                                                                // rho
        static_cast<float>((p_infty/(R*t_infty))*(mach_infty*sqrt(gamma*R*t_infty))),                           // rho*u
        0.0f,                                                                                                   // rho*v
        static_cast<float>((p_infty/(R*t_infty))*(cv*t_infty + 0.5*pow(mach_infty*sqrt(gamma*R*t_infty), 2)))   // rho*E
    };
    for (int i = 0; i < i_max;/*< because nodes->cells*/ i++) {
        q[i][0] = q_inlet;
    }

    for (int j = 1; j < j_max; j++) {
        // inner-lower
            i_bdry = 2;
            i_ghost = 1;
            bdry_velocity = {q[i_bdry][j][1]/q[i_bdry][j][0], q[i_bdry][j][2]/q[i_bdry][j][0]}; // velocity of the boundary cell associated with the inner ghost (outer boundary)
            wall_vec = {mesh_data[2][j+1][0] - mesh_data[2][j][0],
                        mesh_data[2][j+1][1] - mesh_data[2][j][1]}; // vector parallel with the wall element
            wall_norm = {static_cast<float>(-wall_vec[1] / sqrt(pow(wall_vec[0], 2)+pow(wall_vec[1], 2))),
                         static_cast<float>( wall_vec[0] / sqrt(pow(wall_vec[0], 2)+pow(wall_vec[1], 2)))}; // perpendicular, inward-pointing unit vector that is normal to the border wall element
            v_dot_n = bdry_velocity[0]*wall_norm[0] + bdry_velocity[1]*wall_norm[1]; // the magnitude of the velocity component going into the wall       
            q[i_ghost][j] = {
                q[i_bdry][j][0],
                static_cast<float>(q[i_bdry][j][1] - 2*wall_norm[0]*v_dot_n),
                static_cast<float>(q[i_bdry][j][2] - 2*wall_norm[1]*v_dot_n),
                q[i_bdry][j][3]
            };
        // outer-lower
            i_bdry = 3;
            i_ghost = 0;
            bdry_velocity = {q[i_bdry][j][1]/q[i_bdry][j][0], q[i_bdry][j][2]/q[i_bdry][j][0]}; // velocity of the boundary cell associated with the inner ghost (outer boundary)
            v_dot_n = bdry_velocity[0]*wall_norm[0] + bdry_velocity[1]*wall_norm[1]; // the magnitude of the velocity component going into the wall       
            q[i_ghost][j] = {
                q[i_bdry][j][0],
                static_cast<float>(q[i_bdry][j][1] - 2*wall_norm[0]*v_dot_n),
                static_cast<float>(q[i_bdry][j][2] - 2*wall_norm[1]*v_dot_n),
                q[i_bdry][j][3]
            };
        // inner-upper
            i_bdry = i_max-3;
            i_ghost = i_max-2; 
            bdry_velocity = {q[i_bdry][j][1]/q[i_bdry][j][0], q[i_bdry][j][2]/q[i_bdry][j][0]}; // velocity of the boundary cell associated with the inner ghost (outer boundary)
            wall_vec = {mesh_data[i_max-2][j+1][0] - mesh_data[i_max-2][j][0],
                        mesh_data[i_max-2][j+1][1] - mesh_data[i_max-2][j][1]}; // vector parallel with the wall element
            wall_norm = {static_cast<float>( wall_vec[1] / sqrt(pow(wall_vec[0], 2)+pow(wall_vec[1], 2))),
                         static_cast<float>(-wall_vec[0] / sqrt(pow(wall_vec[0], 2)+pow(wall_vec[1], 2)))}; // perpendicular, inward-pointing unit vector that is normal to the border wall element
            v_dot_n = bdry_velocity[0]*wall_norm[0] + bdry_velocity[1]*wall_norm[0];
            q[i_ghost][j] = {
                q[i_bdry][j][0],
                static_cast<float>(q[i_bdry][j][1] - 2*wall_norm[0]*v_dot_n),
                static_cast<float>(q[i_bdry][j][2] - 2*wall_norm[1]*v_dot_n),
                q[i_bdry][j][3]
            };
        // outer-upper
            i_bdry = i_max-4;
            i_ghost = i_max-1; 
            bdry_velocity = {q[i_bdry][j][1]/q[i_bdry][j][0], q[i_bdry][j][2]/q[i_bdry][j][0]}; // velocity of the boundary cell associated with the inner ghost (outer boundary)
            v_dot_n = bdry_velocity[0]*wall_norm[0] + bdry_velocity[1]*wall_norm[0];
            q[i_ghost][j] = {
                q[i_bdry][j][0],
                static_cast<float>(q[i_bdry][j][1] - 2*wall_norm[0]*v_dot_n),
                static_cast<float>(q[i_bdry][j][2] - 2*wall_norm[1]*v_dot_n),
                q[i_bdry][j][3]
            };         
    }


    // outlet boundary condition

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

    boundary_conditions();
    /*
    for (int j = 0; j<j_max; j++) {
        for (int i = 0; i<i_max; i++) {

        }
    }
    */

}
