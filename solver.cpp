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

    float p_infty {static_cast<float>(rho_infty*R*t_infty)};
    float u_infty {static_cast<float>(sqrt(gamma*R*t_infty)*mach_infty)};
    float E_infty {static_cast<float>(0.5*pow(mach_infty * sqrt(gamma*R*t_infty), 2) + (cv*t_infty))};

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
                    static_cast<float>(rho_infty), // solve for density
                    static_cast<float>(rho_infty*u_infty),
                    0.0f,
                    static_cast<float>(rho_infty*E_infty)
                };
                f[i][j] = std::vector<float> {0.0f, 0.0f, 0.0f, 0.0f};
                g[i][j] = std::vector<float> {0.0f, 0.0f, 0.0f, 0.0f};
        }
    }

    std::cout << "Initial values: \nrho: " << rho_infty << "\n" 
    "rho*u: " <<  (rho_infty) * (u_infty)  << "\n" <<
    "rho*V: " << 0.0f << "\n" <<
    "rho*E: " << (rho_infty)*(E_infty) << "\n";
    
    new_q = q;
    update_BCs();
    q = new_q;
}

float Solution::p(int i, int j) {
    // pressure (unused, invalid)
    std::cout << "Pressure getter reporting (shouldnt be called)\n";
    return 0.0;
}     

float Solution::T(int i, int j) {
    // get the energy for the current cell
    float E {static_cast<float>(new_q[i][j][3]/new_q[i][j][0])};
    float V_squared {static_cast<float>(
        pow(new_q[i][j][1]/new_q[i][j][0], 2) + pow(new_q[i][j][2]/new_q[i][j][0], 2)
    )};
    return (E-0.5*V_squared)/cv;
} 

float Solution::rho(int i, int j) {
    // currently unused
    return new_q[i][j][0];
}

float Solution::e(int i, int j) {
    // unused, invalid
    std::cout << "specific internal static energy getter reporting, invalid \n";
    exit(1);
    return 0.0;
}


// internal functions
float Solution::l(int i, int j, float off_i, float off_j) {
    // returns the length of a cell wall given the wall's index in off-integer notation

    // safety check
    if (off_i > 0.2 && off_j > 0.2) {
        std::cout << "Cannot have both i and j be off-integer.\n";
        exit(1);
    }

    // I know it's terrible, but it's easy to diagnose for now.
    if (off_j > 0.1) {
        return sqrt(
        pow(mesh_data[i][j+1][0] - mesh_data[i+1][j+1][0], 2) + 
        pow(mesh_data[i][j+1][1] - mesh_data[i+1][j+1][1], 2));
    }
    else if (off_i > 0.1) {
        return sqrt(
        pow(mesh_data[i+1][j][0] - mesh_data[i+1][j+1][0], 2) + 
        pow(mesh_data[i+1][j][1] - mesh_data[i+1][j+1][1], 2));
    }
    else if (off_j < 0.1) {
        return sqrt(
        pow(mesh_data[i][j][0] - mesh_data[i+1][j][0], 2) + 
        pow(mesh_data[i][j][1] - mesh_data[i+1][j][1], 2));
    }
    else if (off_i < 0.1) {
        return sqrt(
        pow(mesh_data[i][j][0] - mesh_data[i][j+1][0], 2) + 
        pow(mesh_data[i][j][1] - mesh_data[i][j+1][1], 2));
    } else {
        std::cout << "Need to select a wall in order to get its length (both i and j were zero)\n";
        exit(1);
        return 0.0f;
    }

}

float Solution::lambda(float i, float j) {
    /* input: j and i in off-integer notation
     body:
         takes the velocity from the two wall ajacent cells
         takes the normal from the wall's normal

         calculates the average speed of sound between the two cells

     returns:
         lambda value

         key difference between this code and the python - need to calculate lambda based on velocity at the wall (average between cells)
    */
    int cell_i_left = round(i-0.1);  // left and right are used conceptually (could be a/b 1/2, etc) Not literal
    int cell_j_left = round(j-0.1);

    int cell_i_right = round (i+0.1);
    int cell_j_right = round(j+0.1);

    //float a {sqrt(gamma*R*T())}

    // returns the eigenvalue at the cell
    return 0.0;
}

// cizmas style
float Solution::cizmas_lambda(int i, int j, float off_i, float off_j) {
    /* input: j and i in off-integer notation
     body:
         takes the velocity of the original cell between the two cells
         takes l = the wall normal of the velocity
     returns:
         l

         key difference between this code and the python - need to calculate lambda based on velocity at the wall (average between cells)
    */

    float speed_o_sound_sonic {static_cast<float>(sqrt(gamma*R*T(i, j)))};
    std::vector<float> cell_V {static_cast<float>(new_q[i][j][1]/new_q[i][j][0]),
                               static_cast<float>(new_q[i][j][2]/new_q[i][j][0])};

    // I know it's terrible, but it's easy to diagnose for now.
    float dy;
    float dx;

    if (off_j > 0.1) {
        dy = (mesh_data[i+1][j+1][1]-mesh_data[i][j+1][1]);
        dx = (mesh_data[i+1][j+1][0]-mesh_data[i][j+1][0]);
    }
    else if (off_i > 0.1) {
        dy = (mesh_data[i+1][j][1]-mesh_data[i+1][j+1][1]);
        dx = (mesh_data[i+1][j][0]-mesh_data[i+1][j+1][0]);
    }
    else if (off_j < -0.1) {
        dy = (mesh_data[i][j][1]-mesh_data[i+1][j][1]);
        dx = (mesh_data[i][j][0]-mesh_data[i+1][j][0]);
    }
    else /* if (off_i < -0.1)*/ {
        dy = (mesh_data[i][j+1][1]-mesh_data[i][j][1]);
        dx = (mesh_data[i][j+1][0]-mesh_data[i][j][0]);
    }
    std::vector<float> normal = {static_cast<float>(dy/sqrt(dy*dy + dx*dx)),
                                 static_cast<float>(-dx/sqrt(dy*dy + dx*dx))};

    // returns the eigenvalue at the cell
    return static_cast<float>(fabs(cell_V[0]*normal[0] + cell_V[1]*normal[1]) + speed_o_sound_sonic);
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

    // temporary variables for enforcement
    std::vector<float> bdry_velocity(2, 0.0f);
    std::vector<float> wall_vec(2, 0.0f);
    std::vector<float> wall_norm(2, 0.0f);
    float v_dot_n;
    int i_bdry;     // i value of the cell that is mirrored by the ghost cell (boundary cell)
    int i_ghost;    // i value of the ghost cell

    float p_infty {static_cast<float>(rho_infty*R*t_infty)};
    float u_infty {static_cast<float>(sqrt(gamma*R*t_infty)*mach_infty)};
    float E_infty {static_cast<float>(0.5*pow(mach_infty * sqrt(gamma*R*t_infty), 2) + (cv*t_infty))};


    // inlet boundary condition 
    std::vector<float> q_inlet = {
        static_cast<float>(rho_infty),                                                                // rho
        static_cast<float>(rho_infty*u_infty),                           // rho*u
        0.0f,                                                                                                   // rho*v
        static_cast<float>(rho_infty*E_infty)   // rho*E
    };


    for (int i = 0; i < i_max;/*< because nodes->cells*/ i++) {
        for (int j = 0; j<=1; j++) {
            for (int k = 0; k<=3; k++) {
                new_q[i][j][k] = q_inlet[k];
            }
        }
    }

    // no penetration wall boundary condition (using ghost cells)
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

        I am calculating energy based on the last density and velocity, not extrapolated density and velocity
        it looks better this way.
        this is debatable, and jameson agrees. */
    
    float exit_v_mag_squared {0.0f};
    for (int i=2; i<i_max-2; i++) {
        for (int j=j_max-2; j<=j_max-1; j++) {
            for (int k=0; k<=3; k++) {
                new_q[i][j][k] = static_cast<float>(2*new_q[i][j_max-4][k] - new_q[i][j_max-3][k]); 
                // change the j_max-4 for j-2 and j_max-3 for j-1 if you want to make it constant gradient instead of zero gradient
            }
            exit_v_mag_squared = pow(new_q[i][j_max-4][1], 2) + pow(new_q[i][j_max-4][2], 2) / pow(new_q[i][j_max-4][0], 2);
            new_q[i][j][3] = static_cast<float>(p_infty/(gamma-1) + new_q[i][j_max-4][0]*(0.5*exit_v_mag_squared)); // calculate rho*E
        }
    }
}



void Solution::update_f_g(int i, int j) {
    // changed from being based on q to being based on new_q
    // confirmed correct (check 'proofs for fluxes.py')
    float p = (new_q[i][j][3] - 0.5*(new_q[i][j][1]*new_q[i][j][1] + new_q[i][j][2]*new_q[i][j][2])/new_q[i][j][0])*(gamma-1);

    // p calculation debuging statement
    // if (i==2 && j==23) {
    //     std::cout << "p = " << p << "=" << new_q[i][j][3] << "- 0.5*(" << new_q[i][j][1] << "*" << new_q[i][j][1] << "+" << new_q[i][j][2] << "*" << new_q[i][j][2] << ")/" << new_q[i][j][0] << ")*(gamma-1)\n";
    //     std::cout << new_q[i][j][3]*new_q[i][j][1]/new_q[i][j][0] + p*(new_q[i][j][1]/new_q[i][j][0]) << "=" <<
    //     new_q[i][j][3] << "*" << new_q[i][j][1] << "/" << new_q[i][j][0]  << " + " <<  p << "*(" << new_q[i][j][1] << "/" << new_q[i][j][0] << ")\n"; 
    // }

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
    // static cast???

    
}

arrayD3 Solution::get_q() {
    return new_q;
}

void Solution::iterate() {
    /* conduct one iteration*/
    // std::cout << "Iteration " << iteration_count << "\n";
    
    // get iteration alpha constant
    constexpr float alpha_values[] = {0.25f, 0.3333334f, 0.5f, 1.0f};
    float alpha = alpha_values[iteration_count % 4];
    float delta_t = 0.0001f;
    float sum_l_lamb = 0.0f; // $\Sum_0^4 l*\lambda$

    float dx_e;
    float dx_n;
    float dx_w;
    float dx_s;

    float dy_e;
    float dy_n;
    float dy_w;
    float dy_s;

    // allocations for calculation tools
    float area; // unused as of now.
    std::vector<float> res(4, 0);
    std::vector<float> curr_dissipation(4, 0);

    // update all f and g
    for (int i = 0; i<i_max; i++) {
        for (int j = 0; j<j_max; j++) {
            update_f_g(i,j); // calculate all f and g for the iteration (correct)
        }   
    }

    for (int i = 2; i<i_max-2; i++) { // start and end at 2, <i_max-2 because we do not generate a residual for ghost cells (BCs take care of that)
        for (int j = 2; j<j_max-2; j++) {
                        
            // calcualte cell area
            area = static_cast<float>(0.5*((mesh_data[i+1][j+1][0] - mesh_data[i][j][0])    *   (mesh_data[i+1][j][1] - mesh_data[i][j+1][1]) - 
                                    (mesh_data[i+1][j+1][1] - mesh_data[i][j][1])          *   (mesh_data[i+1][j][0] - mesh_data[i][j+1][0])));

            // cell wall deltas (correct and in-line with py script)
            dy_e =  (mesh_data[i+1][j+1][1]-mesh_data[i][j+1][1]);
            dx_e =  (mesh_data[i+1][j+1][0]-mesh_data[i][j+1][0]);

            dy_n =  (mesh_data[i+1][j][1]-mesh_data[i+1][j+1][1]);
            dx_n =  (mesh_data[i+1][j][0]-mesh_data[i+1][j+1][0]);

            dy_w =  (mesh_data[i][j][1]-mesh_data[i+1][j][1]);
            dx_w =  (mesh_data[i][j][0]-mesh_data[i+1][j][0]);

            dy_s =  (mesh_data[i][j+1][1]-mesh_data[i][j][1]);
            dx_s =  (mesh_data[i][j+1][0]-mesh_data[i][j][0]);


            if (j==23 && i==2) {
                std::cout << "------ " << f[i+1][j][0] << "------ \n" <<
                f[i][j-1][0] << "," << f[i][j][0] <<  "," << f[i][j+1][0] << "\n" << 
                "------ " << f[i-1][j][0] << " ------\n";
                // std::cout << "for clarity in the confusion, f[i][j-1][0] was last seen as " << f[i][j-1][0];
                // std::cout << "f[i][j-1][0] = " << f[i][j-1][0] << std::endl;
            }

            // multiply f by dy and g by (-dx). The reason for the x-y and sign reversal is to make the wall delta into a normal
            for (int k = 0; k<=3; k++) {
                res[k] = static_cast<float>(  0.5*(f[i][j][k] + f[i][j+1][k])*dy_e  +  0.5*(g[i][j][k] + g[i][j+1][k])*(-dx_e)        // east
                                             +0.5*(f[i][j][k] + f[i+1][j][k])*dy_n  +  0.5*(g[i][j][k] + g[i+1][j][k])*(-dx_n)        // north
                                             +0.5*(f[i][j][k] + f[i][j-1][k])*dy_w  +  0.5*(g[i][j][k] + g[i][j-1][k])*(-dx_w)        // west
                                             +0.5*(f[i][j][k] + f[i-1][j][k])*dy_s  +  0.5*(g[i][j][k] + g[i-1][j][k])*(-dx_s)        // south
                );

                // // debugging print statement
                // if (j==23 && i==2 && k==2) {
                //     std::cout <<  0.5*(f[i][j][k] + f[i][j+1][k]) << "*" << dy_e << "+" << 0.5*(g[i][j][k] + g[i][j+1][k]) << "*" << (-dx_e) << "\n" << 
                //                  +0.5*(f[i][j][k] + f[i+1][j][k]) << "*" << dy_n << "+" << 0.5*(g[i][j][k] + g[i+1][j][k]) << "*" << (-dx_n) << "\n" << 
                //                  +0.5*(f[i][j][k] + f[i][j-1][k]) << "*" << dy_w << "+" << 0.5*(g[i][j][k] + g[i][j-1][k]) << "*" << (-dx_w) << "\n" << 
                //                  +0.5*(f[i][j][k] + f[i-1][j][k]) << "*" << dy_s << "+" << 0.5*(g[i][j][k] + g[i-1][j][k]) << "*" << (-dx_s) << "\n" <<
                //                  "= " << res[k] << "\n";
                // }
            }

            // calculate dissipation $\vec D$ every 4
            if (iteration_count%4 == 0) {
                curr_dissipation = D(i, j);
            }

            // half index notation is stupid. Henceforth, we will have i,j,i_offset,j_offset
            sum_l_lamb = sqrt(pow(dx_e, 2)+pow(dy_e, 2)) * cizmas_lambda(i, j, 0.0f, 0.5f) +
                         sqrt(pow(dx_n, 2)+pow(dy_n, 2)) * cizmas_lambda(i, j, 0.5f, 0.0f) +
                         sqrt(pow(dx_w, 2)+pow(dy_w, 2)) * cizmas_lambda(i, j, 0.0f, -0.5f) +
                         sqrt(pow(dx_s, 2)+pow(dy_s, 2)) * cizmas_lambda(i, j, -0.5f, 0.0f);

            // new q update debugging statement
            if (i==2 && j==23) {
                //std::cout << "sum_l_lamb: " << sum_l_lamb << "\n";
                for (int k = 0; k<=3; k++) {
                    std::cout << "new q:[" << k << "] = " << q[i][j][k] << "-" << (alpha * CFL * 2 / sum_l_lamb) << "*" << (res[k] - curr_dissipation[k]) << "\n";
                }

            //     std::cout << "lengths\n" <<
            //                 "east " << sqrt(pow(dx_e, 2)+pow(dy_e, 2)) << "\n" <<
            //                 "north" << sqrt(pow(dx_n, 2)+pow(dy_n, 2)) << "\n" <<
            //                 "west" << sqrt(pow(dx_w, 2)+pow(dy_w, 2)) << "\n" <<
            //                 "south" << sqrt(pow(dx_s, 2)+pow(dy_s, 2)) << "\n";
            //     std::cout << "Lambdas\n" <<
            //     "east " << cizmas_lambda(i, j, 0.0f, 0.5f) << "\n" <<
            //     "north " << cizmas_lambda(i, j, 0.5f, 0.0f) << "\n" <<
            //     "west " << cizmas_lambda(i, j, 0.0f, -0.5f)<< "\n" <<
            //     "south " << cizmas_lambda(i, j, -0.5f, 0.0f)     << "\n";
            //     // speed of sound calculation has been confirmed to be the same
            //     // That means either the normal or the velocity is a bit off
            }

            // update the new q
            for (int k = 0; k<4; k++) {
                new_q[i][j][k] = static_cast<float>(q[i][j][k] - (alpha * CFL * 2 / sum_l_lamb) * (res[k] - curr_dissipation[k])); // residual and dissipation
            }
            update_f_g(i,j); // update f and g for the current cell
        }
    }

    update_BCs();           // update the boundary conditions to match inner calculations

    // for (int i = 21; i<32; i++) {
    //     new_q
    // }

    // update q based on intermediate q
    if (iteration_count%4 == 3) {
        q = new_q;              // update inner cells based on inner calculations
    }

    // key difference between this code and the python - is the main q updated every time

    iteration_count++;      // update iteration counter
}
