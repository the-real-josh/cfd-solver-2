#include <iostream>
#include "solver.h"



int main() {
    std::cout << "main reporting\n";

    // get solver configuration (from python)
    get_config();

     // initialize the solver
    arrayD3 mesh_data = get_mesh_data();
    Solution sol;
    sol.innit(mesh_data);

    // run the solver
    for (int i = 0; i<=3; i++) {
        sol.iterate();
    }

    // save results
    save_data(sol.get_q());

    std::cout << "Process finished with exit code 0";
}