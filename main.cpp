#include <iostream>
#include "solver.h"

int main() {
    std::cout << "main reporting";
    get_mesh_data(); // from csv; should come in with ghost cells pre-appended
    initialize(0.3); // based on mesh size
    Solution sol;
    sol.iterate();
    sol.save();
}