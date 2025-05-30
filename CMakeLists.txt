cmake_minimum_required (VERSION 3.30) 
project (CFD_solver_2) # Project name, any name

include(Ctest) # for testing in the future
enable_testing()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# the documentation said to add this to cmakelists.txt
# add_subdirectory(external/csv-parser)

# important - add executable
# specifies executable name (run.exe)
# followed by *all source and header files*
add_executable(run 
            solver.h
            solver.cpp
            setup.cpp
            main.cpp)

# set_property(TARGET CFD_solver_2 PROPERTY CXX_STANDARD 17)

target_include_directories(run PRIVATE external/csv-parser/single_include)
                                 #external\csv-parser\single_include