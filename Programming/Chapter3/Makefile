
EIGEN_FLAX = -I /usr/include/eigen3
# use jsoncpp library
JSONCPP_FLAX = -I /usr/include/jsoncpp

CXX = g++

SRC_DIR = ./src
OBJ_DIR = ./lib
BIN_DIR = ./bin
INCLUDE_DIR = ./include
INCLUDE_DIRS = $(shell find $(INCLUDE_DIR) -type d)
INCLUDE_FLAGS = $(patsubst %,-I%,$(INCLUDE_DIRS))


CXXFLAGS = $(INCLUDE_FLAGS) $(EIGEN_FLAX) $(JSONCPP_FLAX) -std=c++17 -Wall -Wextra -Werror -pedantic -g -flto -fsanitize=address

SRC_FILES = $(shell find $(SRC_DIR) -name '*.cc')

OBJ_FILES = $(patsubst $(SRC_DIR)/%.cc, $(OBJ_DIR)/%.o, $(SRC_FILES))
OBJ_DIRS  = $(dir $(OBJ_FILES))
OBJ_FLAGS = $(patsubst %,-L%,$(OBJ_DIRS))

MAIN_FILES = $(wildcard *.cc)

TARGETS = $(patsubst %.cc, $(BIN_DIR)/%, $(MAIN_FILES))

all: $(OBJ_FILES) $(TARGETS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	@mkdir -p $(OBJ_DIRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGETS): $(BIN_DIR)/%: %.cc $(OBJ_FILES)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< $(OBJ_FILES) -o $@

clean:
	rm -rf $(OBJ_DIR) $(TARGETS) $(BIN_DIR) $(OBJ_DIRS) $(OBJ_FILES)

.PHONY: all clean

