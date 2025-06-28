#include <iostream>
#include "solver.h"
#include <cmath>

int main() {
    // get solver configuration (from python)
    get_config();

    // initialize the solver
    arrayD3 mesh_data = get_mesh_data();
    Solution sol;
    sol.innit(mesh_data);

    // run the solver
    for (int i = 0; i<=max_iterations; i++) {
        sol.iterate();

        // check for convergence every 500
        if (i%500==0 && i>=17000) {
            if (sol.residuals[round(i/4)-1] > (0.95)*sol.residuals[round(i/4)-4000]) {
                // if residuals are not halfing every 150, then it has probably converged
                std::cout<< "Convergence detected. Residual is " << sol.residuals[round(i/4)-1] << "and 800 iterations ago it was " << sol.residuals[round(i/4)-150] << "\n";
                break;
            }
        }
    }

    // save results
    save_data(sol.get_q(), res_fname);
    save_residuals(sol.residuals);
    std::cout << "Process finished with exit code 0\n";
}