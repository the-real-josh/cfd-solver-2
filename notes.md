overall structure:
 - declare functions and classes in the .h file
 - in subordinate files, fill in definitions
 - subordinate files
    - main.cpp - overall structure for the program
    - setup.cpp - get mesh data and initialize with flow
    - solver.cpp - body of the solver

 questions:
 - what links the files
 - modularization
 - structures: equivalent of tuple? has to be declared ahead of time?

makefile is run with make.exe via command prompt. get make.exe with minGW
https://medium.com/@samsorrahman/how-to-run-a-makefile-in-windows-b4d115d7c516

# help

## help on structuring a C++ program
 https://www.cs.odu.edu/~zeil/cs333/f13/Public/cppProgramStructure/cppProgramStructure-web.pdf#page=46.00

## help on writing the solver
https://stackoverflow.com/questions/7910364/beginners-guide-to-own-cfd-code-2d-euler-equation


now we play a game called where is the dissipation error?
Fluxes:
 - f and g calculation (not there, see python script proof)
 - f and g references (unlikely)

Dissipations:
 - switches
 - 

# residuals

Python
alphas[it_number]*2*CFL / sum_l_lamb): 0.25*2*1.0 / 200.1911187852921)
cell lambdas: [('north ', np.float64(369.3528458934361)), ('east ', np.float64(447.52078429176606)), ('south', np.float64(382.3780816400992)), ('west ', np.float64(450.55014898364414))]
new_q[0] = 1.1766 - (0.0024976132959037872)*(-6.893633070859439 - 0.0)
new_q[1] = 122.5619589101473 - (0.0024976132959037872)*(-634.7984179078464 - -19.11058276212613)
new_q[2] = 0.0 - (0.0024976132959037872)*(-229.93326619997586 - 52.76038523344379)
new_q[3] = 259693.21583956483 - (0.0024976132959037872)*(-2115179.0673582414 - 0.0)

x momentum dissipations at i=22,j=2: 0.0 + 0.0 - 0.0 - 19.110582762126143
0.17881806827951885 * 28.43010460657416 + 0.24668861703995185 * -56.86020921314838


C++
Alpha * CFL * 2 / sum_l_lamb = 0.25 * 1 * 2 / 200.347
new_q[2][22][0] = 1.1766 - 0.00249567 * (-6.89362 - 0)
new_q[2][22][1] = 122.562 - 0.00249567 * (-634.796 - -14.679)
new_q[2][22][2] = 0 - 0.00249567 * (-229.933 - 40.5256)
new_q[2][22][3] = 259692 - 0.00249567 * (-2.11517e+06 - 0)

x momentum dissipations at i=22,j=2: 0 + 0 - 0 - 19.1265
28.4301 * 0.178971 + -56.8601 * 0.246893

Analyzing the reason for the similarities in coeffs:
comparing switches: 0.00390625 vs 0.00390625
Comparing lengths 0.123939 vs 0.165157
Comparing lambdas 369.67 vs 382.695

Alpha * CFL * 2 / sum_l_lamb = 0.25 * 1 * 2 / 200.347
new_q[2][22][0] = 1.1766 - 0.00249567 * (-6.89362 - 0)
new_q[2][22][1] = 122.562 - 0.00249567 * (-634.796 - -19.1265)
new_q[2][22][2] = 0 - 0.00249567 * (-229.933 - 52.8045)
new_q[2][22][3] = 259692 - 0.00249567 * (-2.11517e+06 - 0.00106128)


differences:
 - C++ has the same lambda on each side 
 - python 

similarities:
 - same finite difference
