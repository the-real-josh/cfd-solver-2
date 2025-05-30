cmake_minimum_required (VERSION 2.6) # set the version of cmake to 2.6 or better (2.8 was out in 2012 lol)


project (CFD_solver_2) # name


# the documentation said to add this to cmakelists.txt
add_subdirectory(csv-parser)

# important - add executable
# specifies executable name (run.exe)
# followed by *all source and header files*
add_executable(run solver.h solver.cpp setup.cpp main.cpp)


target_link_libraries(run csv)