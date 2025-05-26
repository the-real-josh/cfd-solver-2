#include <iostream> // console messages
#include <vector> // vectors my beloved
#include "external/csv-parser/include/csv.hpp" // csv reading (how to import??)
#include "solver.h"


// setup
arrayD3 get_mesh_data() {
    std::cout << "mesh data gettter reporting\n";

    // get the data
    /*
    csv::CSVReader reader("test.csv");
    for (csv::CSVRow& row : reader) {
        std::string name = row["Name"].get<>();
        int age = row["Age"].get<int>();
    } */
    // return it as an arrayD3
    
    arrayD3 mesh_data;
    return mesh_data;
}


