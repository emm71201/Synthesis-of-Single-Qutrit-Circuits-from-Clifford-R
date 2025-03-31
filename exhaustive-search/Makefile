# # Compiler
# CXX = clang++

# # Compiler Flags
# # CXXFLAGS = -g -O3 -Xpreprocessor -fopenmp -std=c++17 -I/usr/include/boost  
# CXXFLAGS = -g -O3 -Xpreprocessor -fopenmp -std=c++17 -I/usr/include/boost -I/usr/local/include

# # Boost library
# # LIBS = -lboost_system -lboost_filesystem -lomp
# LIBS = -lboost_system -lboost_filesystem -L/usr/local/lib -lomp -L/opt/homebrew/opt/libomp/lib

# # Source Files
# SRCS = main.cpp utils.cpp eiseinstein.cpp cnumber_approximation.cpp 

# # Object Files
# OBJS = $(SRCS:.cpp=.o)

# # Header Files
# DEPS = utils.h eiseinstein.h cnumber_approximation.h 

# # Output executable
# TARGET = process

# # Default rule to build the target
# all: $(TARGET)

# # Rule to link the object files and create the executable
# $(TARGET): $(OBJS)
# 	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

# # Rule to compile each .cpp file into an object file
# %.o: %.cpp $(DEPS)
# 	$(CXX) $(CXXFLAGS) -c $< -o $@

# # Clean up generated files
# clean:
# 	rm -f $(OBJS) $(TARGET)

# # Rebuild the project
# rebuild: clean all

#---------

# # Compiler
# CXX = clang++

# # Compiler Flags
# CXXFLAGS = -g -O3 -Xpreprocessor -fopenmp -std=c++17 -I/usr/include/boost -I/opt/homebrew/opt/libomp/include

# # Boost and OpenMP libraries
# LIBS = -lboost_system -lboost_filesystem -L/usr/local/lib -L/opt/homebrew/opt/libomp/lib -lomp

# # Source Files
# SRCS = main.cpp utils.cpp eiseinstein.cpp cnumber_approximation.cpp 

# # Object Files
# OBJS = $(SRCS:.cpp=.o)

# # Header Files
# DEPS = utils.h eiseinstein.h cnumber_approximation.h 

# # Output executable
# TARGET = process

# # Default rule to build the target
# all: $(TARGET)

# # Rule to link the object files and create the executable
# $(TARGET): $(OBJS)
# 	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

# # Rule to compile each .cpp file into an object file
# %.o: %.cpp $(DEPS)
# 	$(CXX) $(CXXFLAGS) -c $< -o $@

# # Clean up generated files
# clean:
# 	rm -f $(OBJS) $(TARGET)

# # Rebuild the project
# rebuild: clean all

# Compiler
CXX = /opt/homebrew/bin/g++-14

# Compiler Flags
CXXFLAGS = -g -O3 -std=c++17 -I/opt/homebrew/include -I/opt/homebrew/opt/libomp/include -Xpreprocessor -fopenmp

# Boost and OpenMP libraries
LIBS = -lboost_system -lboost_filesystem -L/opt/homebrew/lib -L/opt/homebrew/opt/libomp/lib -lomp

# Source Files
SRCS = main.cpp utils.cpp eiseinstein.cpp cnumber_approximation.cpp 

# Object Files
OBJS = $(SRCS:.cpp=.o)

# Header Files
DEPS = utils.h eiseinstein.h cnumber_approximation.h 

# Output executable
TARGET = process

# Default rule to build the target
all: $(TARGET)

# Rule to link the object files and create the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

# Rule to compile each .cpp file into an object file
%.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up generated files
clean:
	rm -f $(OBJS) $(TARGET)

# Rebuild the project
rebuild: clean all


