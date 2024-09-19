CXX = g++

SRC_DIR = ./src
OBJ_DIR = ./lib
BIN_DIR = .
INCLUDE_DIR = ./include

CXXFLAGS = -I$(INCLUDE_DIR) -std=c++11 -Wall -O2

SRC_FILES = $(wildcard $(SRC_DIR)/*.cc)

OBJ_FILES = $(patsubst $(SRC_DIR)/%.cc, $(OBJ_DIR)/%.o, $(SRC_FILES))

MAIN_FILES = $(wildcard *.cc)

TARGETS = $(patsubst %.cc, $(BIN_DIR)/%, $(MAIN_FILES))

all: $(OBJ_FILES) $(TARGETS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGETS): $(BIN_DIR)/%: %.cc $(OBJ_FILES)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< $(OBJ_FILES) -o $@

clean:
	rm -f $(TARGETS)

distclean: clean
	rm -f $(OBJ_DIR)/*

run: $(TARGETS)
	@for target in $(TARGETS); do \
		echo "Running $$target";  \
		./$$target;               \
		echo "\n";                \
	done

# Phony targets
.PHONY: all clean distclean run
