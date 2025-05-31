@echo off
echo Compiling project...

g++ -std=c++17 -Wall -Wextra -O2 -c solver.cpp
g++ -std=c++17 -Wall -Wextra -O2 -c setup.cpp
g++ -std=c++17 -Wall -Wextra -O2 -c main.cpp
::do not need to compile csv.hpp because the library is a header-only library??

g++ -o main.exe main.o setup.o solver.o

echo Done.
