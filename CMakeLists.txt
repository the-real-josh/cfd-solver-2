cmake_minimum_required (VERSION 3.30) 
project (solver_engine) # Project name, any name

include(Ctest) # for testing in the future
enable_testing()

# todo: consistent naming

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# the documentation said to add this to cmakelists.txt
# add_subdirectory(external/csv-parser)

# important - add executable
# specifies executable name (run.exe)
# followed by *all source and header files*
add_executable(solver_engine
            solver.h
            solver.cpp
            setup.cpp
            main.cpp)

# set_property(TARGET solver_engine PROPERTY CXX_STANDARD 17)

target_include_directories(solver_engine PRIVATE external/csv-parser/single_include)
                                 #external\csv-parser\single_include