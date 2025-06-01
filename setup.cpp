#include <iostream> // console messages
#include <vector> // vectors my beloved
#include "csv.hpp" // do not need paths because that is taken care of in the includepaths in cmake
#include "solver.h"

// official declaration of config globals goes here.
int i_max;
int j_max;
float t_infty; // kelvin
float p_infty; // pascals
float mach_infty;
std::string mesh_fname;
std::string res_fname;

// setup
void get_config() {
    // get these guys based on a config file later
    i_max = 3;
    j_max = 3;
    t_infty = 300.0f;
    p_infty = 101325.0f;
    mach_infty = 0.3f;
    mesh_fname = "mesh sh=11x21.csv";
    res_fname = "output.csv";
}

arrayD3 get_mesh_data() {
    std::cout << "mesh data gettter reporting\n";

    // TODO: find out the dimensions from config file or cin
    // np.size-like. ie, an array {1, 1, 1} has n_cols of 3 and n_rows of 0

    // get the data
    arrayD2 flat_mesh_data;

    // read csv in flat-paired form
    csv::CSVReader reader("test.csv");
    for (auto& row: reader) {
        // Note: Can also use index of column with [] operator
        flat_mesh_data.push_back({row["x"].get<float>(), row["y"].get<float>()});
    }

    // // debug print
    // for (const auto& pair: flat_mesh_data) {
    //     for (const auto& xy: pair) {
    //         std::cout << xy << " "; 
    //     }
    //     std::cout << " -- ";
    // }
    // std::cout << "\n";
    
    
    // reshape into proper dimensions

    // declare mesh_data
    arrayD3 mesh_data;

    // assert that the given dimensions match the data
    std::cout << "mesh data size: " << flat_mesh_data.size() << "\n";
    std::cout << "mesh integer division: " << flat_mesh_data.size() % i_max << "\n";

    //if(i_max == 0 || flat_mesh_data.size()%i_max != 0 ) throw std::domain_error( "bad #cols" ) ;

    // fill up mesh_data
    for (int j = 0; j < j_max; j++) {
        // one whole row each iteration using 
        mesh_data.push_back(
            arrayD2(flat_mesh_data.begin() + (j*i_max),
                    flat_mesh_data.begin() + ((j+1)*i_max))
        );
    }

    // // print out the mesh_data to see if it is good
    // for (const auto& row : mesh_data) {
    //     for (const auto& pair: row) {
    //         for (const auto& xy: pair) {
    //             std::cout << xy;
    //         }
    //         std::cout << " ";
    //     }
    //     std::cout << "\n";
    // }

    // return it as an arrayD3
    return mesh_data;
}


void save_data(arrayD3 q_out) {
    /* save the data
       saving format: 4 columns for each of the state variables */

    // flatten the cells
    std::vector<float> flat_rho;
    std::vector<float> flat_rho_u;
    std::vector<float> flat_rho_v;
    std::vector<float> flat_rho_E;
    
    for (int j = 0; j<q_out.size(); j++) {
        for (int i = 0; i< q_out[0].size(); i++) {
            flat_rho.push_back(q_out[j][i][0]);
            flat_rho_u.push_back(q_out[j][i][1]);
            flat_rho_v.push_back(q_out[j][i][2]);
            flat_rho_E.push_back(q_out[j][i][3]);
        }
    }
    
    
    // write to the output file
    std::ofstream file(res_fname);
    auto writer = csv::make_csv_writer(file);
    writer << std::vector<std::string>({"rho", "rho_u", "rho_v", "rho_E"});
    for (int i = 0; i<flat_rho.size(); i++) {
        writer << std::vector<std::string>({std::to_string(flat_rho[i]),
             std::to_string(flat_rho_u[i]),
             std::to_string(flat_rho_v[i]),
             std::to_string(flat_rho_E[i])
            });  
    }

}

