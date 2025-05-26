# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

# Target executable
TARGET = main

# Source files
SRCS = main.cpp setup.cpp solver.cpp

# Object files (replace .cpp with .o)
OBJS = $(SRCS:.cpp=.o)

# Default target
all: $(TARGET)

# Link object files into the final executable
$(TARGET): $(OBJS)
$(CXX) $(CXXFLAGS) -o $@ $^

# Compile .cpp files into .o files
%.o: %.cpp solver.h
$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build artifacts
clean:
rm -f $(OBJS) $(TARGET)
