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
    iteration_count = 0;

    std::cout << "i_max, j_max = " << i_max << "," << j_max << "\n";
    // system("pause");

    /*initialize arrays
        q gets the inlet boundary condition duplicated throughout the domain (uniform rightward flow)
        f and g get zeros of correct shape*/
    q.resize(i_max);
    f.resize(i_max);
    g.resize(i_max);
    for (int i = 0; i<=i_max-1; i++) { 
        q[i].resize(j_max);
        f[i].resize(j_max);
        g[i].resize(j_max);
        for (int j = 0; j<=j_max-1; j++) {
                q[i][j] = std::vector<float>{
                    static_cast<float>(p_infty / (R * t_infty)), // solve for density
                    static_cast<float>((p_infty / (R * t_infty)) * (mach_infty * sqrt(gamma*R*t_infty))),
                    0.0f,
                    static_cast<float>((p_infty/(R*t_infty))*(0.5*pow(mach_infty * sqrt(gamma*R*t_infty), 2) + (cv*t_infty)))
                };
                f[i][j] = std::vector<float> {0.0f, 0.0f, 0.0f, 0.0f};
                g[i][j] = std::vector<float> {0.0f, 0.0f, 0.0f, 0.0f};
        }
    }
    update_BCs();
    new_q = q;
}

// simple state getters
float Solution::p(int i, int j) {
    std::cout << "Pressure getter reporting";
    return 0.0;
}                 
float Solution::T(int i, int j) {
    std::cout << "temperature getter reporting";
    return 0.0;
}                 
float Solution::rho(int i, int j) {
    std::cout << "density getter reporting";
    return 0.0;
}
float Solution::e(int i, int j) {
    std::cout << "specific internal static energy getter reporting";
    return 0.0;
}


// internal functions
float Solution::l(float i, float j) {
    // returns the length of a cell wall given the wall's index in off-integer notation
    std::cout << "length function reporting";
    // convert to wall-coordinates
    return 0.0;
    // returns the length of a cell wall  
}

float Solution::lambda(float i, float j) {
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

std::vector<float> Solution::D(int i, int j) {
    float nu_2 {0.25};
    float nu_4 {0.002};
    /*input: 
        j and i cells
    body;
        calculates second order dissipations
        calculates fourth order dissipations
    returns:
        total dissipations*/
    
    return std::vector<float> {0.0f, 0.0f, 0.0f, 0.0f};
}

void Solution::update_BCs() {
    /*input:
        None
    body: 
        edits q such that boundary conditions are enforced
    output: 
        None*/
    std::cout << "Enforcing Boundary Conditions\n";

    // temporary variables for enforcement
    std::vector<float> bdry_velocity(2, 0.0f);
    std::vector<float> wall_vec(2, 0.0f);
    std::vector<float> wall_norm(2, 0.0f);
    float v_dot_n;
    int i_bdry;     // i value of the cell that is mirrored by the ghost cell (boundary cell)
    int i_ghost;    // i value of the ghost cell


    // inlet boundary condition 
    std::vector<float> q_inlet = {
        static_cast<float>(p_infty/(R*t_infty)),                                                                // rho
        static_cast<float>((p_infty/(R*t_infty))*(mach_infty*sqrt(gamma*R*t_infty))),                           // rho*u
        0.0f,                                                                                                   // rho*v
        static_cast<float>((p_infty/(R*t_infty))*(cv*t_infty + 0.5*pow(mach_infty*sqrt(gamma*R*t_infty), 2)))   // rho*E
    };

    for (int i = 0; i < i_max;/*< because nodes->cells*/ i++) {
        for (int j = 0; j<=1; j++) {
            for (int k = 0; k<=3; k++) {
                q[i][j][k] = q_inlet[k];
            }
        }
        // needs to be j=0 and j=1, not just j=0
    }

    // no penetration wall boundary condition (using ghost cells)
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
    /*Energy boundary condition:
        do not apply zero gradient to energy
        If subsonic:
            get the exit pressure boundary condition (idk make it same as inlet)
            If supersonic: Calculate it with the pressure right before the end of the domain (i_max - )
            
        the python script made both ghost cells have the same value
        
        In cizmas's notes, he said that q'' is zero, not q', so it would imply that they have a constant "slope". Ask the TA.

        In the Jameson paper, he said to say p = p_infty, and to do something with eigenvalues? He also says that there are other outer boundary conditions 
        that have been proposed that are worth investigating. I guess if it isn't playing tennis with the inlet boundary then it's fine? 
        */
    // no not apply zero gradient to energy. Calculate it with the exit pressure (subsonic) or the 
    float exit_v_mag_squared {0.0f};
    for (int i=2; i<i_max-2; i++) {
        for (int j=i_max-2; j<=i_max-1; j++) {
            for (int k=0; k<=3; k++) {
                q[i][j][k] = static_cast<float>(2*q[i][j_max-4][k] - q[i][j_max-3][k]); 
                // change the j_max-4 for j-2 and j_max-3 for j-1 if you want to make it constant gradient instead of zero gradient
            }
            exit_v_mag_squared = pow(q[i][j][1], 2) + pow(q[i][j][2], 2) / pow(q[i][j][0], 2);
            q[i][j][3] = static_cast<float>(p_infty / (gamma-1) + q[i][j][0]*(0.5*exit_v_mag_squared)); // calculate rho*E
        }
    }
}

void Solution::update_f() {
    for (int i; i<i_max; i++) {
        for (int j; j<j_max; j++) {
            f[i][j][0] = static_cast<float>(q[i][j][1]);
            f[i][j][1] = static_cast<float>(q[i][j][1]/q[i][j][0] + q[i][j][3]*(gamma-1));
            f[i][j][2] = static_cast<float>(pow(q[i][j][1], 2)*q[i][j][2]/q[i][j][0]);
            f[i][j][3] = static_cast<float>(q[i][j][3]*q[i][j][1]*gamma/q[i][j][0]);
        }
    }
}

void Solution::update_g() {
    for (int i; i<i_max; i++) {
        for (int j; j<j_max; j++) {
            g[i][j][0] = static_cast<float>(q[i][j][2]);
            g[i][j][1] = static_cast<float>(q[i][j][1]*q[i][j][2]/q[i][j][0]);
            g[i][j][2] = static_cast<float>(pow(q[i][j][2], 2)/q[i][j][0] + q[i][j][3]*(gamma-1));
            g[i][j][3] = static_cast<float>(q[i][j][3]*q[i][j][2]*gamma / q[i][j][0]);
        }
    }
}

arrayD3 Solution::get_q() {
    return q;
}

void Solution::iterate() {
    /* conduct one iteration*/
    std::cout << "Iteration " << iteration_count << "\n";
    
    // get iteration alpha constant
    constexpr float alpha_values[] = {0.25f, 0.3333334f, 0.5f, 1.0f};
    float alpha = alpha_values[iteration_count % 4];
    float delta_t = 0.0001f;

    // allocations for calculation tools
    float area;
    std::vector<float> res(4, 0);
    std::vector<float> curr_dissipation(4, 0);


    // calculate all f and g for the iteration
    update_f(); 
    update_g();

    for (int i = 2; i<i_max-2; i++) { // start and end at 2, <i_max-2 because we do not generate a residual for ghost cells (BCs take care of that)
        for (int j = 2; j<j_max-2; j++) {
        
        std::cout << "\n\nDEBUG: begin analysis of cell i,j= " << i << "," << j << "\n";
        
        // calcualte cell area
        area = static_cast<float>(0.5*((mesh_data[i+1][j+1][0] - mesh_data[i][j][0])    *   (mesh_data[i+1][j][1] - mesh_data[i][j+1][1]) - 
                                 (mesh_data[i+1][j+1][1] - mesh_data[i][j][1])          *   (mesh_data[i+1][j][0] - mesh_data[i][j+1][0])));
        // std::cout << "calculating area based on the cells with coords \n" << 
        // mesh_data[i+1][j+1][0] << "," << mesh_data[i+1][j+1][1] << "\n" << 
        // mesh_data[i+1][j][0] << "," << mesh_data[i+1][j][1] << "\n" << 
        // mesh_data[i][j][0] << "," << mesh_data[i][j][1] << "\n" << 
        // mesh_data[i][j+1][0] << "," << mesh_data[i][j+1][1] << "\n" <<
        // "resulting area: " << area << "\n";

        // calculate residual
        // residual = fs*delta ys, gs*delta xs
        for (int k = 0; k<3; k++) {
            res[k] = 
                static_cast<float>(0.5*(f[i][j][k] + f[i+1][j][k])*(mesh_data[i+1][j+1][1]-mesh_data[i][j+1][1])  -   0.5*(g[i][j][k] + g[i+1][j][k])*(mesh_data[i+1][j+1][0]-mesh_data[i][j+1][0])) +       // i+ in f and g, and i+1/2 in dy and dx
                static_cast<float>(0.5*(f[i][j][k] + f[i][j+1][k])*(mesh_data[i+1][j+1][1]-mesh_data[i+1][j][1])  -   0.5*(g[i][j][k] + g[i][j+1][k])*(mesh_data[i+1][j+1][0]-mesh_data[i+1][j][0])) +       // j+
                static_cast<float>(0.5*(f[i][j][k] + f[i-1][j][k])*(mesh_data[i+1][j][1]-mesh_data[i][j][1])      -   0.5*(g[i][j][k] + g[i-1][j][k])*(mesh_data[i+1][j][0]-mesh_data[i][j][0])) +           // i-
                static_cast<float>(0.5*(f[i][j][k] + f[i][j-1][k])*(mesh_data[i+1][j][1]-mesh_data[i][j][1])      -   0.5*(g[i][j][k] + g[i][j+1][k])*(mesh_data[i+1][j+1][0]-mesh_data[i+1][j][0]));        // j-
        }

        // calculate dissipation $\vec D$
        curr_dissipation = D(i, j);
        for (auto &p : curr_dissipation) {
            std::cout << p;
        }

        // update the new q
        for (int k = 0; k<3; k++) {
            new_q[i][j][k] = CFL * delta_t * (1/area) * (res[k] - curr_dissipation[k]);
        }

        // calculate time_step
        // update new_q
        // system("pause");
        }
    }

    q = new_q;              // update inner cells based on inner calculations
    update_BCs();           // update the boundary conditions to match inner calculations
    iteration_count++;      // update iteration counter
}
