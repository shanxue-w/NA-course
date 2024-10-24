CXX = g++

SRC_DIR = ./src
OBJ_DIR = ./lib
BIN_DIR = .
INCLUDE_DIR = ./include
INCLUDE_DIRS = $(shell find $(INCLUDE_DIR) -type d)
INCLUDE_FLAGS = $(patsubst %,-I%,$(INCLUDE_DIRS))

CXXFLAGS = $(INCLUDE_FLAGS) -std=c++11 -Wall -O3

SRC_FILES = $(shell find $(SRC_DIR) -name '*.cc')

OBJ_FILES = $(patsubst $(SRC_DIR)/%.cc, $(OBJ_DIR)/%.o, $(SRC_FILES))
OBJ_DIRS  = $(dir $(OBJ_FILES))

MAIN_FILES = $(wildcard *.cc)

TARGETS = $(patsubst %.cc, $(BIN_DIR)/%, $(MAIN_FILES))

all: $(OBJ_FILES) $(TARGETS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	@mkdir -p $(OBJ_DIRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGETS): $(BIN_DIR)/%: %.cc $(OBJ_FILES)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< $(OBJ_FILES) -o $@

report: HW2_3220104819_Wanghao.tex 
	# check the file ./figures/B_Newton2.png exists, then xelatex
	# if not, run the program to generate the file
	if [ ! -f ./figures/B_Newton2.png ]; then \
		make run && make plot; \
	fi
	xelatex HW2_3220104819_Wanghao.tex 
	xelatex HW2_3220104819_Wanghao.tex 

run: $(TARGETS)
	@for target in $(TARGETS); do \
		echo "Running $$target";  \
		./$$target;               \
		echo "\n";                \
	done

./data/B_Newton2.txt: $(TARGETS)
	make run
plot: ./data/B_Newton2.txt
	python draw.py

doxygen: 
	doxygen Doxyfile
	cd latex && make

clean:
	rm -f $(TARGETS) *.log *.aux *.bbl *.blg *.out *.pdf *.run.xml *.toc *.bib
	rm -rf $(OBJ_DIR)/* latex/ html/ data/* 

distclean: clean
	rm -rf figures/*

# Phony targets
.PHONY: all clean distclean run report doxygen

