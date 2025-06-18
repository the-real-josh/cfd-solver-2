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

 