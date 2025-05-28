#include <iostream> // console messages
#include <vector> // vectors my beloved
#include "external/csv-parser/include/csv.hpp" // csv reading (how to import??)
#include "solver.h"

// setup
arrayD3 create_3d_array(int depth, int rows, int cols) {
    return arrayD3(depth, std::vector<std::vector<float>>(rows, std::vector<float>(cols, 0.0f)));
}

// something is bad here.
arrayD3 get_mesh_data() {
    std::cout << "mesh data gettter reporting\n";

    // TODO: find out the dimensions from commands.txt
    //int dimensions[2] = {3, 3};

    // get the data
    arrayD2 flat_mesh_data(2);
    //arrayD3 mesh_data = create_3d_array(dimensions[0], dimensions[1], 4);

    // read csv
    csv::CSVReader reader("test.csv");
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
    // reshape into proper dimensions


    // return it as an arrayD3
    arrayD3 mesh_data;
    return mesh_data;
}


// save the data
void save_data(arrayD3) {
    // SOMEHOW save the data
}

