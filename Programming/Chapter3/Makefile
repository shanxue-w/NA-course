
EIGEN_FLAX = -I/usr/include/eigen3
# use jsoncpp library
JSONCPP_FLAX = -I/usr/include/jsoncpp -I/usr/include/jsoncpp/json -L/usr/local/lib -ljsoncpp
LAPACK_FLAX = -DEIGEN_USE_LAPACKE -llapacke -llapack -lblas -I/usr/include/lapacke
GMP_FLAX = -lgmpxx -lgmp -L/usr/local/lib -I/usr/local/include

CXX = g++

SRC_DIR = ./src
OBJ_DIR = ./lib
BIN_DIR = ./bin
INCLUDE_DIR = ./include
INCLUDE_DIRS = $(shell find $(INCLUDE_DIR) -type d)
INCLUDE_FLAGS = $(patsubst %,-I%,$(INCLUDE_DIRS))
ADDRESS_FLAG = -fsanitize=address

CXXFLAGS = $(INCLUDE_FLAGS) $(GMP_FLAX) $(EIGEN_FLAX) $(JSONCPP_FLAX) -std=c++17 -Wall -Wextra -Werror -pedantic -O3 -fuse-linker-plugin -mavx -mfma -march=native -fPIC -funroll-loops -fprefetch-loop-arrays -fomit-frame-pointer -fopenmp 

SRC_FILES = $(shell find $(SRC_DIR) -name '*.cc')

OBJ_FILES = $(patsubst $(SRC_DIR)/%.cc, $(OBJ_DIR)/%.o, $(SRC_FILES))
OBJ_DIRS  = $(dir $(OBJ_FILES))
OBJ_FLAGS = $(patsubst %,-L%,$(OBJ_DIRS))

MAIN_FILES = $(wildcard *.cc)

TARGETS = $(patsubst %.cc, $(BIN_DIR)/%, $(MAIN_FILES))

all: $(OBJ_FILES) $(TARGETS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	@mkdir -p $(OBJ_DIRS)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

$(TARGETS): $(BIN_DIR)/%: %.cc $(OBJ_FILES)
	@mkdir -p $(BIN_DIR)
	$(CXX) $< $(OBJ_FILES) -o $@ $(CXXFLAGS)

doxygen:
	doxygen Doxyfile
	cd latex && make

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR) latex/ html/


.PHONY: all clean

