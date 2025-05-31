cmake_minimum_required (VERSION 3.30) 
project (solver_engine) # Project name, any name

include(Ctest) # for testing in the future
enable_testing()


set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# important - add executable
# specifies executable name (run.exe)
# followed by *all source and header files*
add_executable(solver_engine
            solver.h
            solver.cpp
            setup.cpp
            main.cpp)

target_include_directories(solver_engine PRIVATE external/csv-parser/single_include)
                                 #external\csv-parser\single_include