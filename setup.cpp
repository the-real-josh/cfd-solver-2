#include <iostream> // console messages
#include <vector> // vectors my beloved
#include "csv.hpp"
#include "solver.h"


// setup
arrayD3 create_3d_array(int depth, int rows, int cols) {
    return arrayD3(depth, std::vector<std::vector<float>>(rows, std::vector<float>(cols, 0.0f)));
}

// something is bad here.
arrayD3 get_mesh_data() {
    std::cout << "mesh data gettter reporting\n";

    // TODO: find out the dimensions from config file or cin
    // np.size-like. ie, an array {1, 1, 1} has n_cols of 3 and n_rows of 0
    int i_max {3}; 
    int j_max {3};

    // get the data
    arrayD2 flat_mesh_data(2);
    //arrayD3 mesh_data = create_3d_array(dimensions[0], dimensions[1], 4);

    // read csv
   /* csv::CSVReader reader("test.csv");
    for (auto& row: reader) {
        // Note: Can also use index of column with [] operator
        flat_mesh_data[0].push_back(row["x"].get<float>());
        flat_mesh_data[1].push_back(row["y"].get<float>());
    }

    for (const auto& i : flat_mesh_data[0]) {
        std::cout << i << " ";
    }
    for (const auto& i : flat_mesh_data[1]) {
        std::cout << i << " ";
    }
    std::cout << "\n";
    
    trying to rework dimensionality here*/

    // read csv
    csv::CSVReader reader("test.csv");
    for (auto& row: reader) {
        // Note: Can also use index of column with [] operator
        flat_mesh_data.push_back({row["x"].get<float>()});
    }

    for (const auto& pair: flat_mesh_data) {
        for (const auto& xy: pair) {
            std::cout << xy << " "; 
        }
    }
    std::cout << "\n";
    
    
    
    
    // reshape into proper dimensions

    // declare mesh_data
    arrayD3 mesh_data;

    // assert that the given dimensions match the data
    std::cout << "mesh data size" << flat_mesh_data.size() << "\n";
    std::cout << "mesh integer division" << flat_mesh_data.size() % i_max << "\n";

    //if(i_max == 0 || flat_mesh_data.size()%i_max != 0 ) throw std::domain_error( "bad #cols" ) ;

    // fill up mesh_data
    for (int j = 0; j < j_max; j++) {
        // one whole row each iteration using 
        mesh_data.push_back(
            arrayD2(flat_mesh_data.begin() + (j*i_max),
                    flat_mesh_data.begin() + ((j+1)*i_max))
        );
    }


    // print out the mesh_data to see if it is good
    for (const auto& row : mesh_data) {
        for (const auto& pair: row) {
            for (const auto& xy: pair) {
                std::cout << xy;
            }
            std::cout << " ";
        }
        std::cout << "\n";
    }



    // return it as an arrayD3
    return mesh_data;

    
}


// save the data
void save_data(arrayD3) {
    // SOMEHOW save the data
}

