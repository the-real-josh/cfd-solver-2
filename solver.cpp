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
                    static_cast<float>((p_infty / (R * t_infty)) * (mach_infty * sqrt(1.4*R*t_infty))),
                    0.0f,
                    static_cast<float>(0.5*pow(mach_infty * sqrt(1.4*R*t_infty), 2) + (cv*t_infty))
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
    std::vector<float> boundary_velocity;
    std::vector<float> wall_vec;
    std::vector<float> wall_norm;
    float v_dot_n;

    // inlet boundary condition (compiles and runs good)
    std::vector<float> q_inlet = {
        static_cast<float>(p_infty/(R*t_infty)),
        static_cast<float>(mach_infty*sqrt(gamma*R*t_infty)),
        0.0f,
        static_cast<float>(cv*t_infty + 0.5*pow(mach_infty*sqrt(gamma*R*t_infty), 2))
    };
    for (int i = 0; i < i_max;/*< because nodes->cells*/ i++) {
        q[i][0] = q_inlet;
    }

    // wall boundary condition (inner ghost cells)
    for (int j = 1; j < j_max; j++) {
        // lower
            boundary_velocity = {q[1][j][1]/q[1][j][0], q[1][j][2]/q[1][j][0]}; // velocity of the boundary cell associated with the inner ghost (outer boundary)
            
            std::cout << "\nwall boundary velocity at 1," << j << "\n"; 
            for (auto& el: boundary_velocity) {
                std::cout << el << "\n";
            }


            wall_vec = {mesh_data[2][j+1][0] - mesh_data[2][j][0],
                                             mesh_data[2][j+1][1] - mesh_data[2][j][1]}; // vector parallel with the wall element

            std::cout << "\nwall element tangent vector at 2," << j+1 << "\n"; 
            for (auto& el: wall_vec) {
                std::cout << el << "\n";
            }

            wall_norm = {static_cast<float>(-wall_vec[1] / sqrt(pow(wall_vec[0], 2)+pow(wall_vec[1], 2))),
                                            static_cast<float>( wall_vec[0] / sqrt(pow(wall_vec[0], 2)+pow(wall_vec[1], 2)))}; // perpendicular, inward-pointing unit vector that is normal to the border wall element
            v_dot_n = boundary_velocity[0]*wall_norm[0] + boundary_velocity[1]*wall_norm[1];
            
            std::cout << "\nwall normal vector at 2," << j << "\n"; 
            for (auto& el: wall_norm) {
                std::cout << el << "\n";
            }
            
            q[1][j] = {
                q[2][j][0],
                static_cast<float>(q[2][j][1] + 2*wall_norm[0]*v_dot_n),
                static_cast<float>(q[2][j][2] + 2*wall_norm[1]*v_dot_n),
                q[2][j][3]
            };

            system("pause");

        // upper
            boundary_velocity = {q[i_max-3][j][1]/q[i_max-3][j][0], q[i_max-3][j][2]/q[i_max-3][j][0]}; // velocity of the boundary cell associated with the inner ghost (outer boundary)
            wall_vec = {mesh_data[i_max-2][j+1][0] - mesh_data[i_max-2][j][0],
                                             mesh_data[i_max-2][j+1][1] - mesh_data[i_max-2][j][1]}; // vector parallel with the wall element
            wall_norm = {static_cast<float>( wall_vec[1] / sqrt(pow(wall_vec[0], 2)+pow(wall_vec[1], 2))),
                                            static_cast<float>(-wall_vec[0] / sqrt(pow(wall_vec[0], 2)+pow(wall_vec[1], 2)))}; // perpendicular, inward-pointing unit vector that is normal to the border wall element
            v_dot_n = boundary_velocity[0]*wall_norm[0] + boundary_velocity[1]*wall_norm[0];
            q[i_max-2][j] = {
                q[i_max-3][j][0],
                static_cast<float>(q[i_max-3][j][1] + 2*wall_norm[0]*v_dot_n),
                static_cast<float>(q[i_max-3][j][2] + 2*wall_norm[1]*v_dot_n),
                q[i_max-3][j][3]
            };

    }

    // wall boundary condition (outer ghost cells)

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
    
    std::cout << "input i, then enter, then j then enter to get mesh data.";
    int test_j;
    int test_i;
    for (int i = 0; i<10; i++) {
        std::cin >> test_i;
        std::cin >> test_j;
        std::cout << "mesh_data for those indeces: ";
        for (auto& k : mesh_data[test_i][test_j]) {
            std::cout << k << "   ";
        }
        std::cout <<  "\n";
    }

    boundary_conditions();
    /*
    for (int j = 0; j<j_max; j++) {
        for (int i = 0; i<i_max; i++) {

        }
    }
    */

}
