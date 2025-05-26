#include <iostream> // console messages
#include <vector> // vectors my beloved
#include "csv.hpp" // csv reading (how to import??)

#define arrayD3 std::vector<std::vector<std::vector<float>>> // float array of dimension 3

// define all the mesh stuff here (empty)

arrayD3 get_mesh_data() {
    // get the data
    csv::CSVReader reader("data.csv");
    for (csv::CSVRow& row : reader) {
        std::string name = row["Name"].get<>();
        int age = row["Age"].get<int>();
    }
    // return it as an arrayD3
}

void allocate() {
    // pre-allocate
}

void run() {
    // run the solver
}


int main() {
    std::cout << "hello cruel world";

    get_mesh_data(); // from csv
    allocate(); // based on mesh size
    run();
}
