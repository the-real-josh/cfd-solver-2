#include <iostream> 
#include <vector>   
#include "csv.hpp"  
#include "solver.h" 
#include <cmath>

void Solution::innit(arrayD3 _mesh_data) {
    std::cout << "innit\n";

    // take the mesh data into the class
    mesh_data = _mesh_data;
    iteration_count = 0;

    // output the mesh size
    std::cout << "i_max, j_max = " << i_max << "," << j_max << "\n";

    // calculate other freestream variables from the provided
    float p_infty {static_cast<float>(rho_infty*R*t_infty)};
    float u_infty {static_cast<float>(sqrt(gamma*R*t_infty)*mach_infty)};
    float E_infty {static_cast<float>(0.5*pow(mach_infty * sqrt(gamma*R*t_infty), 2) + (cv*t_infty))};

    /*initialize arrays:
        q gets uniform rightward flow
        f, g fluxes and D dissipation get zeros of correct shape*/
    q.resize(i_max);
    f.resize(i_max);
    g.resize(i_max);
    D.resize(i_max);
    for (int i = 0; i<=i_max-1; i++) { 
        q[i].resize(j_max);
        f[i].resize(j_max);
        g[i].resize(j_max);
        D[i].resize(j_max);
        for (int j = 0; j<=j_max-1; j++) {
                q[i][j] = std::vector<float>{
                    static_cast<float>(rho_infty), // solve for density
                    static_cast<float>(rho_infty*u_infty),
                    0.0f,
                    static_cast<float>(rho_infty*E_infty)
                };
                f[i][j] = std::vector<float> {0.0f, 0.0f, 0.0f, 0.0f};
                g[i][j] = std::vector<float> {0.0f, 0.0f, 0.0f, 0.0f};
                D[i][j] = std::vector<float> {0.0f, 0.0f, 0.0f, 0.0f};
        }
    }

    // print freestream values
    std::cout << "Initial values: \nrho: " << rho_infty << "\n" 
    "rho*u: " <<  (rho_infty) * (u_infty)  << "\n" <<
    "rho*V: " << 0.0f << "\n" <<
    "rho*E: " << (rho_infty)*(E_infty) << "\n";
    
    // ensure that the boundary conditions are enforced in initialization for both q and new_q
    new_q = q;
    update_BCs();
    q = new_q;
}

float Solution::p(int i, int j) {
    // calculate pressure at cell i,j
    // Used in the dissipation calculation
    return static_cast<float>((new_q[i][j][3] - 0.5*(new_q[i][j][1]*new_q[i][j][1] + new_q[i][j][2]*new_q[i][j][2])/new_q[i][j][0])*(gamma-1));
}     

float Solution::T(int i, int j) {
    // Returns the temperature for the cell i,j
    float E {static_cast<float>(new_q[i][j][3]/new_q[i][j][0])};
    float V_squared {static_cast<float>(
        pow(new_q[i][j][1]/new_q[i][j][0], 2) + pow(new_q[i][j][2]/new_q[i][j][0], 2)
    )};
    return (E-0.5*V_squared)/cv;
} 


// internal functions
float Solution::l(int i, int j, float off_i, float off_j) {
    /* returns the length of a cell wall given the wall's index in off-integer index notation*/

    float dy;
    float dx;

    if (fabs(off_j-0.5) < small && fabs(off_i) < small) { // 
        dy = (mesh_data[i+1][j+1][1]-mesh_data[i][j+1][1]);
        dx = (mesh_data[i+1][j+1][0]-mesh_data[i][j+1][0]);
    }
    else if (fabs(off_i-0.5) < small && fabs(off_j) < small) {
        dy = (mesh_data[i+1][j][1]-mesh_data[i+1][j+1][1]);
        dx = (mesh_data[i+1][j][0]-mesh_data[i+1][j+1][0]);
    }
    else if (fabs(off_j+0.5) < small && fabs(off_i) < small) {
        dy = (mesh_data[i][j][1]-mesh_data[i+1][j][1]);
        dx = (mesh_data[i][j][0]-mesh_data[i+1][j][0]);
    }
    else if (fabs(off_i+0.5) < small && fabs(off_j) < small) {
        dy = (mesh_data[i][j+1][1]-mesh_data[i][j][1]);
        dx = (mesh_data[i][j+1][0]-mesh_data[i][j][0]);
    } else {
        std::cout << "Please select a wall to get the eigenvalue at.";
        exit(1);
    }

    return static_cast<float>(sqrt(dy*dy + dx*dx));
}

float Solution::lambda(int i, int j, float off_i, float off_j) {
    /*Returns the fastest wave speed at i,j in a direction specified by off_i, off_j*/
    
    float speed_of_sound {static_cast<float>(sqrt(gamma*R*T(i, j)))};

    std::vector<float> cell_V {static_cast<float>(new_q[i][j][1]/new_q[i][j][0]),
                               static_cast<float>(new_q[i][j][2]/new_q[i][j][0])};

    float dy;
    float dx;

    if (fabs(off_j-0.5) < small && fabs(off_i) < small) { // 
        dy = (mesh_data[i+1][j+1][1]-mesh_data[i][j+1][1]);
        dx = (mesh_data[i+1][j+1][0]-mesh_data[i][j+1][0]);
    }
    else if (fabs(off_i-0.5) < small && fabs(off_j) < small) {
        dy = (mesh_data[i+1][j][1]-mesh_data[i+1][j+1][1]);
        dx = (mesh_data[i+1][j][0]-mesh_data[i+1][j+1][0]);
    }
    else if (fabs(off_j+0.5) < small && fabs(off_i) < small) {
        dy = (mesh_data[i][j][1]-mesh_data[i+1][j][1]);
        dx = (mesh_data[i][j][0]-mesh_data[i+1][j][0]);
    }
    else if (fabs(off_i+0.5) < small && fabs(off_j) < small) {
        dy = (mesh_data[i][j+1][1]-mesh_data[i][j][1]);
        dx = (mesh_data[i][j+1][0]-mesh_data[i][j][0]);
    } else {
        std::cout << "Please select a wall to get the eigenvalue at.";
        exit(1);
    }

    std::vector<float> normal = {static_cast<float>(dy/sqrt(dy*dy + dx*dx)),
                                 static_cast<float>(-dx/sqrt(dy*dy + dx*dx))};

    // returns the eigenvalue at the cell - abs(Dot product of velocity with cell wall normal), plus speed of sound
    return static_cast<float>(fabs(cell_V[0]*normal[0] + cell_V[1]*normal[1]) + speed_of_sound);
}


float Solution::switch_2_xi(int i, int j, float off_i, float off_j) {
    /*Second order switch - detect sharp fluctuations in pressure, then activate the
    second order dissipation.

    The switch is requested at the wall but calculated at the cell. Therefore it is calculated
    at the two wall-bordering cells, then averaged.*/
    if (fabs(off_i) > 0.01) {
        std::cout << "Error in switch_2_xi";
        exit(1);
    }

    // original cell
    float central_switch = nu_2 * fabs(p(i,j+1) - 2*p(i,j) + p(i,j-1)) / 
           (p(i,j+1) + 2*p(i,j) + p(i,j-1));
    
    // bordering cell (confirmed to find a different cell every time)
    int border_i = round(i+2*off_i);
    int border_j = round(j+2*off_j);

    float border_switch = nu_2 * fabs(p(border_i,border_j+1) - 2*p(border_i,border_j) + p(border_i,border_j-1)) / 
           (p(border_i,border_j+1) + 2*p(border_i,border_j) + p(border_i,border_j-1));

    return 0.5f*(central_switch + border_switch); // resulting value is always positve
}
float Solution::switch_2_eta(int i, int j, float off_i, float off_j) {
    /*Second order switch - detect sharp fluctuations in pressure, then activate the
    second order dissipation.

    The switch is requested at the wall but calculated at the cell. Therefore it is calculated
    at the two wall-bordering cells, then averaged.
    
    The variation in i in the numerator and variation in j in the denominator is NOT a typo.*/

    if (fabs(off_j) > 0.01) {
        std::cout << "Error in switch_2_eta";
        exit(1);
    }

    // original cell
    float central_switch =  nu_2 * fabs(p(i+1, j) - 2*p(i, j) + p(i-1, j)) / 
            (p(i,j+1) + 2*p(i,j) + p(i,j-1));

    // bordering cell (confirmed to find a different cell every time)
    int border_i = round(i+2*off_i);
    int border_j = round(j+2*off_j); 

    float border_switch =  nu_2 * fabs(p(border_i+1, border_j) - 2*p(border_i, border_j) + p(border_i-1, border_j)) / 
        (p(border_i,border_j+1) + 2*p(border_i,border_j) + p(border_i,border_j-1));

    return 0.5f*(central_switch + border_switch); // resulting value is always postive
}


float Solution::switch_4_xi(int i, int j, float off_i, float off_j) {
    /*4th order switch - return nu_4, unless the 2nd order switch is large.
    Then, return 0 to disable 4th order dissipation*/ 
    float sw = nu_4 - switch_2_xi(i, j, off_i, off_j);
    if (sw < 0.0f) {
        return 0.0f;
    } else {
        return sw;
    }
}
float Solution::switch_4_eta(int i, int j, float off_i, float off_j) {
    /*4th order switch - return nu_4, unless the 2nd order switch is large.
    Then, return 0 to disable 4th order dissipation*/ 
    float sw = nu_4 - switch_2_eta(i, j, off_i, off_j);
    if (sw < 0.0f) {
        return 0.0f;
    } else {
        return sw;
    }
}

std::vector<float> Solution::update_D(int i, int j) {
    /* Gives the correct dissipation value at any given cell i,j
        Applies Jameson second and fourth order dissipation terms to damp out the euler-instability wiggles and aid convergence on a solution.*/

    std::vector<float> final_dissipation(4, 0.0f);

    // 2nd order coefficients
    float coeff_1 = switch_2_xi(i,j,   0.0f,  0.5f) * l(i,j,  0.0f,  0.5f) * lambda(i,j,  0.0f,  0.5f);
    float coeff_2 = switch_2_xi(i,j,   0.0f, -0.5f) * l(i,j,  0.0f, -0.5f) * lambda(i,j,  0.0f, -0.5f);

    float coeff_3 = switch_2_eta(i,j,  0.5f,  0.0f) * l(i,j,  0.5f,  0.0f) * lambda(i,j,  0.5f,  0.0f);
    float coeff_4 = switch_2_eta(i,j, -0.5f,  0.0f) * l(i,j, -0.5f,  0.0f) * lambda(i,j, -0.5f,  0.0f);


    // 4th order coefficients
    float coeff_5 = switch_4_xi(i,j,  0.0f,  0.5f) * l(i,j,  0.0f,  0.5f) * lambda(i,j,  0.0f,  0.5f);
    float coeff_6 = switch_4_xi(i,j,  0.0f, -0.5f) * l(i,j,  0.0f, -0.5f) * lambda(i,j,  0.0f, -0.5f);

    float coeff_7 = switch_4_eta(i,j, 0.5f,  0.0f) * l(i,j,  0.5f,  0.0f) * lambda(i,j,  0.5f,  0.0f);
    float coeff_8 = switch_4_eta(i,j,-0.5f,  0.0f) * l(i,j, -0.5f,  0.0f) * lambda(i,j, -0.5f,  0.0f);

    // multiply coefficients by finite differences
    for (int k = 0; k<4; k++) {
        final_dissipation[k] = static_cast<float>(
                                    // 2nd order dissipation terms
                                    (new_q[i][j+1][k] - new_q[i][j+0][k])*coeff_1  -
                                    (new_q[i][j+0][k] - new_q[i][j-1][k])*coeff_2  +
 
                                    (new_q[i+1][j][k] - new_q[i+0][j][k])*coeff_3  -
                                    (new_q[i+0][j][k] - new_q[i-1][j][k])*coeff_4  +

                                    // 4th order dissipation terms
                                    -(new_q[i][j+2][k] - 3*new_q[i][j+1][k] + 3*new_q[i][j+0][k] - new_q[i][j-1][k])*coeff_5 - 
                                    -(new_q[i][j+1][k] - 3*new_q[i][j+0][k] + 3*new_q[i][j-1][k] - new_q[i][j-2][k])*coeff_6

                                    -(new_q[i+2][j][k] - 3*new_q[i+1][j][k] + 3*new_q[i+0][j][k] - new_q[i-1][j][k])*coeff_7 - 
                                    -(new_q[i+1][j][k] - 3*new_q[i+0][j][k] + 3*new_q[i-1][j][k] - new_q[i-2][j][k])*coeff_8
                                );
    }
    return final_dissipation;
}

void Solution::update_BCs() {
    /* edits new_q such that boundary conditions are enforced */

    // temporary variables
    std::vector<float> bdry_velocity(2, 0.0f);  // velocity of real cell ajacent with wall
    std::vector<float> wall_vec(2, 0.0f);       // vector that is tangent with wall border
    std::vector<float> wall_norm(2, 0.0f);      // vector that is normal to the wall border
    float v_dot_n;  // component of the velocity that is normal to the wall
    int i_bdry;     // i value of the cell that is mirrored by the ghost cell (boundary cell)
    int i_ghost;    // i value of the ghost cell

    // freestream variables for reference
    float p_infty {static_cast<float>(rho_infty*R*t_infty)};
    float u_infty {static_cast<float>(sqrt(gamma*R*t_infty)*mach_infty)};
    float E_infty {static_cast<float>(0.5*pow(mach_infty * sqrt(gamma*R*t_infty), 2) + (cv*t_infty))};

    // inlet boundary condition 
    std::vector<float> q_inlet = {
        static_cast<float>(rho_infty),          // rho
        static_cast<float>(rho_infty*u_infty),  // rho*u
        0.0f,                                   // rho*v
        static_cast<float>(rho_infty*E_infty)   // rho*E
    };
    for (int i = 0; i < i_max;/*< because nodes->cells*/ i++) {
        for (int j = 0; j<=1; j++) {
            for (int k = 0; k<=3; k++) {
                new_q[i][j][k] = q_inlet[k];
            }
        }
    }


    // outlet boundary condition
    float exit_v_mag_squared {0.0f};
    for (int i=2; i<i_max-2; i++) {
        for (int j=j_max-2; j<=j_max-1; j++) {
            for (int k=0; k<=3; k++) {
                // calculate energy boundary condition based on constant gradient
                new_q[i][j][k] = static_cast<float>(2*new_q[i][j-1][k] - new_q[i][j-2][k]); 
            }
            // calculate energy separately: based on exit pressure and current cell values of q1 and q2
            exit_v_mag_squared = (pow(new_q[i][j][1], 2) + pow(new_q[i][j][2], 2)) / pow(new_q[i][j][0], 2);
            new_q[i][j][3] = static_cast<float>(p_infty/(gamma-1) + new_q[i][j][0]*(0.5f*exit_v_mag_squared)); // calculate rho*E
        }
    }


    // No penetration wall boundary condition (using ghost cells)
    for (int j = 1; j < j_max; j++) {
        // inner-lower
            i_bdry = 2;
            i_ghost = 1;
            bdry_velocity = {new_q[i_bdry][j][1]/new_q[i_bdry][j][0], new_q[i_bdry][j][2]/new_q[i_bdry][j][0]}; // velocity of the boundary cell associated with the inner ghost (outer boundary)
            wall_vec = {mesh_data[2][j+1][0] - mesh_data[2][j][0],
                        mesh_data[2][j+1][1] - mesh_data[2][j][1]}; // vector parallel with the wall element
            wall_norm = {static_cast<float>(-wall_vec[1] / sqrt(pow(wall_vec[0], 2)+pow(wall_vec[1], 2))),
                         static_cast<float>( wall_vec[0] / sqrt(pow(wall_vec[0], 2)+pow(wall_vec[1], 2)))}; // perpendicular, inward-pointing unit vector that is normal to the border wall element
            v_dot_n = bdry_velocity[0]*wall_norm[0] + bdry_velocity[1]*wall_norm[1]; // the magnitude of the velocity component going into the wall       
            new_q[i_ghost][j] = {
                new_q[i_bdry][j][0],
                static_cast<float>(new_q[i_bdry][j][0]*(new_q[i_bdry][j][1]/new_q[i_bdry][j][0] - 2*wall_norm[0]*v_dot_n)),
                static_cast<float>(new_q[i_bdry][j][0]*(new_q[i_bdry][j][2]/new_q[i_bdry][j][0] - 2*wall_norm[1]*v_dot_n)),
                new_q[i_bdry][j][3]
            };
        // outer-lower
            i_bdry = 3;
            i_ghost = 0;
            bdry_velocity = {new_q[i_bdry][j][1]/new_q[i_bdry][j][0], new_q[i_bdry][j][2]/new_q[i_bdry][j][0]}; // velocity of the boundary cell associated with the inner ghost (outer boundary)
            v_dot_n = bdry_velocity[0]*wall_norm[0] + bdry_velocity[1]*wall_norm[1]; // the magnitude of the velocity component going into the wall       
            new_q[i_ghost][j] = {
                new_q[i_bdry][j][0],
                static_cast<float>(new_q[i_bdry][j][0]*(new_q[i_bdry][j][1]/new_q[i_bdry][j][0] - 2*wall_norm[0]*v_dot_n)),
                static_cast<float>(new_q[i_bdry][j][0]*(new_q[i_bdry][j][2]/new_q[i_bdry][j][0] - 2*wall_norm[1]*v_dot_n)),
                new_q[i_bdry][j][3]
            };
        // inner-upper
            i_bdry = i_max-3;
            i_ghost = i_max-2; 
            bdry_velocity = {new_q[i_bdry][j][1]/new_q[i_bdry][j][0], new_q[i_bdry][j][2]/new_q[i_bdry][j][0]}; // velocity of the boundary cell associated with the inner ghost (outer boundary)
            wall_vec = {mesh_data[i_max-2][j+1][0] - mesh_data[i_max-2][j][0],
                        mesh_data[i_max-2][j+1][1] - mesh_data[i_max-2][j][1]}; // vector parallel with the wall element
            wall_norm = {static_cast<float>( wall_vec[1] / sqrt(pow(wall_vec[0], 2)+pow(wall_vec[1], 2))),
                         static_cast<float>(-wall_vec[0] / sqrt(pow(wall_vec[0], 2)+pow(wall_vec[1], 2)))}; // perpendicular, inward-pointing unit vector that is normal to the border wall element
            v_dot_n = bdry_velocity[0]*wall_norm[0] + bdry_velocity[1]*wall_norm[1];
            new_q[i_ghost][j] = {
                new_q[i_bdry][j][0],
                static_cast<float>(new_q[i_bdry][j][0]*(new_q[i_bdry][j][1]/new_q[i_bdry][j][0] - 2*wall_norm[0]*v_dot_n)),
                static_cast<float>(new_q[i_bdry][j][0]*(new_q[i_bdry][j][2]/new_q[i_bdry][j][0] - 2*wall_norm[1]*v_dot_n)),
                new_q[i_bdry][j][3]
            };
        // outer-upper
            i_bdry = i_max-4;
            i_ghost = i_max-1; 
            bdry_velocity = {new_q[i_bdry][j][1]/new_q[i_bdry][j][0], new_q[i_bdry][j][2]/new_q[i_bdry][j][0]}; // velocity of the boundary cell associated with the inner ghost (outer boundary)
            v_dot_n = bdry_velocity[0]*wall_norm[0] + bdry_velocity[1]*wall_norm[1];
            new_q[i_ghost][j] = {
                new_q[i_bdry][j][0],
                static_cast<float>(new_q[i_bdry][j][0]*(new_q[i_bdry][j][1]/new_q[i_bdry][j][0] - 2*wall_norm[0]*v_dot_n)),
                static_cast<float>(new_q[i_bdry][j][0]*(new_q[i_bdry][j][2]/new_q[i_bdry][j][0] - 2*wall_norm[1]*v_dot_n)),
                new_q[i_bdry][j][3]
            };        
    }
}

void Solution::update_f_g(int i, int j) {
    /*Updates the vectors f and g based on new_q*/
    
    float p = (new_q[i][j][3] - 0.5*(new_q[i][j][1]*new_q[i][j][1] + new_q[i][j][2]*new_q[i][j][2])/new_q[i][j][0])*(gamma-1);

    // update f
    f[i][j][0] = static_cast<float>(new_q[i][j][1]);
    f[i][j][1] = static_cast<float>(new_q[i][j][1]*new_q[i][j][1]/new_q[i][j][0] + p);
    f[i][j][2] = static_cast<float>(new_q[i][j][1]*new_q[i][j][2]/new_q[i][j][0]);
    f[i][j][3] = static_cast<float>(new_q[i][j][3]*new_q[i][j][1]/new_q[i][j][0] + p*(new_q[i][j][1]/new_q[i][j][0]));

    // update g
    g[i][j][0] = static_cast<float>(new_q[i][j][2]);
    g[i][j][1] = static_cast<float>(new_q[i][j][1]*new_q[i][j][2]/new_q[i][j][0]);
    g[i][j][2] = static_cast<float>(new_q[i][j][2]*new_q[i][j][2]/new_q[i][j][0] + p); 
    g[i][j][3] = static_cast<float>(new_q[i][j][3]*new_q[i][j][2]/new_q[i][j][0] + p*(new_q[i][j][2]/new_q[i][j][0])); // for clarity in the confusion, f[i][j-1][0] was last seen as 122.5620*0.0775929+101324*0.0158403 ??????
}

arrayD3 Solution::get_q() {
    /* read-only access of current state for external functions */
    return new_q;
}

void Solution::iterate() {
    /* conduct one iteration*/
    
    // get iteration alpha constant
    constexpr float alpha_values[] = {0.25f, 0.3333334f, 0.5f, 1.0f};
    float alpha = alpha_values[iteration_count % 4];
    float sum_l_lamb = 0.0f; // $\Sum_0^4 l*\lambda$

    float ag_res = 0.0f; // aggregate residual

    // for cell normals
    float dx_e;
    float dx_n;
    float dx_w;
    float dx_s;
    float dy_e;
    float dy_n;
    float dy_w;
    float dy_s;

    // Cell residual
    std::vector<float> res(4, 0.0f);

    // update all f and g
    for (int i = 0; i<i_max; i++) {
        for (int j = 0; j<j_max; j++) {
            update_f_g(i,j); // calculate all f and g for the iteration (correct)
        }   
    }

    // Print iteration count every 100
    if (iteration_count%100 == 0) {
        std::cout << "Iteration " << iteration_count << "\n";
    }

    // residual calculation for all cells
    for (int i = 2; i<i_max-2; i++) { // start and end at 2, <i_max-2 because we do not generate a residual for ghost cells (BCs take care of that)
        for (int j = 2; j<j_max-2; j++) {
                        
            // cell wall deltas: East, North, West, South, counterclockwise
            dy_e =  (mesh_data[i+1][j+1][1]-mesh_data[i][j+1][1]);
            dx_e =  (mesh_data[i+1][j+1][0]-mesh_data[i][j+1][0]);

            dy_n =  (mesh_data[i+1][j][1]-mesh_data[i+1][j+1][1]);
            dx_n =  (mesh_data[i+1][j][0]-mesh_data[i+1][j+1][0]);

            dy_w =  (mesh_data[i][j][1]-mesh_data[i+1][j][1]);
            dx_w =  (mesh_data[i][j][0]-mesh_data[i+1][j][0]);

            dy_s =  (mesh_data[i][j+1][1]-mesh_data[i][j][1]);
            dx_s =  (mesh_data[i][j+1][0]-mesh_data[i][j][0]);


            // multiply f by dy and g by (-dx). The reason for the x-y and sign reversal is to make the wall delta into a normal
            for (int k = 0; k<=3; k++) {
                res[k] = static_cast<float>(  0.5*(f[i][j][k] + f[i][j+1][k])*dy_e  +  0.5*(g[i][j][k] + g[i][j+1][k])*(-dx_e)        // east
                                             +0.5*(f[i][j][k] + f[i+1][j][k])*dy_n  +  0.5*(g[i][j][k] + g[i+1][j][k])*(-dx_n)        // north
                                             +0.5*(f[i][j][k] + f[i][j-1][k])*dy_w  +  0.5*(g[i][j][k] + g[i][j-1][k])*(-dx_w)        // west
                                             +0.5*(f[i][j][k] + f[i-1][j][k])*dy_s  +  0.5*(g[i][j][k] + g[i-1][j][k])*(-dx_s)        // south
                                            );
                if (res[k] != res[k]) {
                    save_data(new_q, res_fname);
                    std::cout << "Nan entry detected at " << i << ", " << j << ", " << k << "\n";
                    exit(1);
                } 
            }

            // RE-calculate dissipation D every 4 iterations. 
            // curr_dissipation is a class-wide 3D array, so it will persist between iterations.
            if (iteration_count%4 == 0) {
                D[i][j] = update_D(i, j);
            }

            // Value for calculation of time-step
            sum_l_lamb = sqrt(pow(dx_e, 2)+pow(dy_e, 2)) * lambda(i, j, 0.0f, 0.5f) +
                         sqrt(pow(dx_n, 2)+pow(dy_n, 2)) * lambda(i, j, 0.5f, 0.0f) +
                         sqrt(pow(dx_w, 2)+pow(dy_w, 2)) * lambda(i, j, 0.0f, -0.5f) +
                         sqrt(pow(dx_s, 2)+pow(dy_s, 2)) * lambda(i, j, -0.5f, 0.0f);

            // update the new_q state vector
            for (int k = 0; k<4; k++) {
                new_q[i][j][k] = static_cast<float>(q[i][j][k] - (alpha * CFL * 2 / sum_l_lamb) * (res[k] - D[i][j][k])); // residual and dissipation
            }
        }
    }

    update_BCs();

    // After a full iteration (4 minor iterations)
    if (iteration_count%4 == 3) {
        // residual calculation
        for (int i = 0; i<i_max-2; i++) {
            for (int j = 0; j<j_max-2; j++) {
                for (int k = 0; k<=3; k++) {
                    ag_res += fabs(q[i][j][k] - new_q[i][j][k]);
                }
            }
        }
        ag_res /= (i_max-2)*(j_max-2);
        residuals.push_back(ag_res);

        // Update official state vector q, after one full iteration (4 minor iterations)
        q = new_q;             
    }

    iteration_count++;      // update iteration counter
}
