@echo off
echo Compiling project...

g++ -std=c++17 -Wall -Wextra -O2 -c solver.cpp
g++ -std=c++17 -Wall -Wextra -O2 -c setup.cpp
g++ -std=c++17 -Wall -Wextra -O2 -c main.cpp

g++ -o main.exe main.o setup.o solver.o

echo Done. Run with: main.exe

pause
