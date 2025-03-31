# Compiler and flags
CXX = /opt/homebrew/bin/g++-14
CXXFLAGS = -std=c++11 -O3 -Iheaders -I/opt/homebrew/Cellar/gmp/6.3.0/include -I/opt/homebrew/Cellar/mpfr/4.2.1/include 
LDFLAGS = -L/opt/homebrew/Cellar/gmp/6.3.0/lib -L/opt/homebrew/Cellar/mpfr/4.2.1/lib -lgmp -lgmpxx -lmpfr

# Target executable
TARGET = main

# Source files
SRC := src/main.cpp src/eiseinstein.cpp src/norm_eq_solver.cpp src/state.cpp src/approximate_cnumber.cpp src/enumeration.cpp  src/state_approximation.cpp
HEADERS := headers/eiseinstein.h src/norm_eq_solver.h headers/state.h headers/approximate_cnumber.h headers/enumeration.h headers/state_approximation.h

# Object files
# The object files must be in the objects directory. If the directory does not exist, it should be created.
OBJ_DIR = objects
OBJ = $(addprefix $(OBJ_DIR)/, $(notdir $(SRC:.cpp=.o)))

# Default rule to build the target
all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $(TARGET) $(LDFLAGS)

# Rule to compile source files into the objects directory
$(OBJ_DIR)/%.o: src/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule for rebuilding (clean + all)
rebuild: clean all

# Clean up generated files
clean:
	rm -rf $(OBJ_DIR) $(TARGET)

# Phony targets
.PHONY: all clean rebuild

