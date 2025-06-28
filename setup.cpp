#include <iostream> // console messages
#include <vector> 
#include "csv.hpp" // do not need paths because that is taken care of in the includepaths in cmake
#include "solver.h"

// official declaration of config globals goes here.
int i_max;
int j_max;
int max_iterations;
float t_infty;      // kelvin
float rho_infty;    // pascals
float mach_infty;
std::string mesh_fname;
std::string res_fname;

// setup
void get_config() {
    // Extract values from a config file
    csv::CSVReader reader("solver_commands.csv"); // hard-coded to be the same filename as the one python makes
    auto it = reader.begin();
    csv::CSVRow row = *it;

    i_max = row["i_max"].get<int>();
    j_max = row["j_max"].get<int>();
    max_iterations = row["max_iterations"].get<int>();
    t_infty = row["t_infty"].get<float>();
    rho_infty = row["rho_infty"].get<float>();
    mach_infty = row["mach_infty"].get<float>();
    mesh_fname = row["mesh_fname"].get<std::string>(); 
    res_fname = row["res_fname"].get<std::string>();


    std::cout << "Config:\n";
    std::cout << "i_max: " << i_max << "\n";
    std::cout << "j_max: " << j_max << "\n";
    std::cout << "Mesh fname: " << mesh_fname << "\n";
}

arrayD3 get_mesh_data() {
    std::cout << "Retrieving Mesh Data\n";

    arrayD2 flat_mesh_data;
    arrayD3 mesh_data;

    // read csv in flat-paired form
    csv::CSVReader reader(mesh_fname);
    for (auto& row: reader) {
        // Note: Can also use index of column with [] operator
        flat_mesh_data.push_back({row["x"].get<float>(), row["y"].get<float>()});
    }

    // reshape into proper dimensions

    // assert that the given dimensions match the data
    std::cout << "mesh data size: " << flat_mesh_data.size() << "\n";
    std::cout << "mesh integer division: " << flat_mesh_data.size() % (i_max+1) << "\n";
    if(flat_mesh_data.size() % (i_max+1) != 0 ) throw std::domain_error( "bad #cols" ) ;

    /* why <= and i_max+1?   see readme.txt#about cell indexing##converting. You are going from max index of nodes to max number of nodes. */
    for (int i = 0; i <= i_max; i++) {
        // one whole row each iteration using 
        mesh_data.push_back(
            arrayD2(flat_mesh_data.begin() + (i*(j_max+1)),
                    flat_mesh_data.begin() + ((i+1)*(j_max+1)))
        );
    }

    return mesh_data;
}

void save_data(arrayD3 q_out, std::string out_filename) {
    /* save the data
       saving format: 4 columns for each of the state variables */

    // flatten the cells
    std::vector<float> flat_rho;
    std::vector<float> flat_rho_u;
    std::vector<float> flat_rho_v;
    std::vector<float> flat_rho_E;
    
    for (int i = 0; i<q_out.size(); i++) {
        for (int j = 0; j< q_out[0].size(); j++) {
            flat_rho.push_back(q_out[i][j][0]);
            flat_rho_u.push_back(q_out[i][j][1]);
            flat_rho_v.push_back(q_out[i][j][2]);
            flat_rho_E.push_back(q_out[i][j][3]);
        }
    }

    // Status message
    std::cout << "Data saver: data has dimensions of i,j " << q_out.size() << "," << q_out[0].size() << "\n";
    
    // write to the output file
    std::ofstream file(out_filename);
    auto writer = csv::make_csv_writer(file);
    writer << std::vector<std::string>({"rho", "rho_u", "rho_v", "rho_E"});
    for (int i = 0; i<flat_rho.size(); i++) {
        writer << std::vector<std::string>({std::to_string(flat_rho[i]),
             std::to_string(flat_rho_u[i]),
             std::to_string(flat_rho_v[i]),
             std::to_string(flat_rho_E[i])
            });  
    }

    std::cout << "Data saver: data has been saved\n";
}


void save_residuals(std::vector<float> residuals) {
    std::ofstream file ("residuals.csv");
    auto writer = csv::make_csv_writer(file);
    writer << std::vector<std::string>({"residuals"});
    for (int i = 0; i<residuals.size(); i++) {
        writer << std::vector<std::string>({std::to_string(residuals[i])});  
    }
}

