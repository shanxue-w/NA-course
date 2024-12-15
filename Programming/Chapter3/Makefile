
EIGEN_FLAX = -I/usr/include/eigen3
# use jsoncpp library
JSONCPP_FLAX = -I/usr/include/jsoncpp -I/usr/include/jsoncpp/json -L/usr/local/lib -ljsoncpp
LAPACK_FLAX = -DEIGEN_USE_LAPACKE -llapacke -llapack -lblas -I/usr/include/lapacke
GMP_FLAX = -lgmpxx -lgmp -L/usr/local/lib -I/usr/local/include

CXX = g++

SRC_DIR = ./src/src
OBJ_DIR = ./lib
BIN_DIR = ./bin
INCLUDE_DIR = ./src/include
INCLUDE_DIRS = $(shell find $(INCLUDE_DIR) -type d)
INCLUDE_FLAGS = $(patsubst %,-I%,$(INCLUDE_DIRS))
ADDRESS_FLAG = -fsanitize=address

CXXFLAGS = $(INCLUDE_FLAGS) $(GMP_FLAX) $(EIGEN_FLAX) $(JSONCPP_FLAX) -std=c++17 -Wall -Wextra -Werror -pedantic -O2 -fuse-linker-plugin -mavx -mfma -fPIC -fopenmp 

SRC_FILES = $(shell find $(SRC_DIR) -name '*.cc')

OBJ_FILES = $(patsubst $(SRC_DIR)/%.cc, $(OBJ_DIR)/%.o, $(SRC_FILES))
OBJ_DIRS  = $(dir $(OBJ_FILES))
OBJ_FLAGS = $(patsubst %,-L%,$(OBJ_DIRS))

# MAIN_FILES = $(wildcard *.cc)
MAIN_DIR = ./src
MAIN_FILES = $(wildcard $(MAIN_DIR)/*.cc)
TARGETS = $(patsubst $(MAIN_DIR)/%.cc, $(BIN_DIR)/%, $(MAIN_FILES))
#echo $(TARGETS)
# MAIN_FILES = $(shell find ./src -name "*.cc")

# TARGETS = $(patsubst ./src/%.cc, $(BIN_DIR)/%, $(MAIN_FILES))

all: $(OBJ_FILES) $(TARGETS) run plot doxygen

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	@mkdir -p $(OBJ_DIRS)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

# $(TARGETS): $(BIN_DIR)/%: %.cc $(OBJ_FILES)
# 	@echo $(TARGETS)
# 	@mkdir -p $(BIN_DIR)
# 	$(CXX) $< $(OBJ_FILES) -o $@ $(CXXFLAGS)

$(TARGETS): $(BIN_DIR)/%: $(MAIN_DIR)/%.cc $(OBJ_FILES)
	@mkdir -p $(BIN_DIR)
	$(CXX) $< $(OBJ_FILES) -o $@ $(CXXFLAGS)

doxygen: plot
	doxygen Doxyfile
	cd doc/design_pdf && make && mv refman.pdf design.pdf && make clean
	cd doc && make && make clean

run: $(TARGETS)
	@mkdir -p result; 
	for target in $(TARGETS); do \
        $$target; \
    done

plot: $(TARGETS) run
	@mkdir -p figure
	@python3 ./src/plot.py

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

.PHONY: all clean

