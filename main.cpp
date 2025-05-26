#include <iostream>
#include "solver.h"


// define inlet state in public namespace
float t_infty = 300; // kelvin
float p_infty = 101325; // pascals


int main() {

    std::cout << "main reporting\n";
     // from csv; should come in with ghost cells pre-appended
    arrayD3 mesh_data = get_mesh_data();
    Solution sol;
    sol.innit(mesh_data, 0.3);
    sol.iterate();
    sol.get_q();
}