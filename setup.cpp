#include <iostream> // console messages
#include <vector> // vectors my beloved
#include "csv.hpp" // csv reading (how to import??)
#include "solver.h"


// setup
arrayD3 get_mesh_data() {

    std::cout << "mesh data gettter reporting";

    // get the data
    /*
    csv::CSVReader reader("data.csv");
    for (csv::CSVRow& row : reader) {
        std::string name = row["Name"].get<>();
        int age = row["Age"].get<int>();
    } */
    // return it as an arrayD3
    arrayD3 mesh_data;
    return mesh_data;
}

void boundary_conditions(float mach) {

    std::cout << "boundary condition getter reporting";
    
    // define inlet state
    float t_infty = 300; // kelvin
    float p_infty = 101325; // pascals
    
    // inlet BC
    // outlet BC
    // wall BC
}

void initialize(float mach) {

    std::cout << "initializer reporting";

    // pre-allocate, depending on the size
    boundary_conditions(mach);
    // initialize cells with all columns identical
}


