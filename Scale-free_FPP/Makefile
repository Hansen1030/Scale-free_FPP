CXX := clang++ 
CXXFLAGS := -Wall -std=c++11 

SRC_DIR := .
BUILD_DIR := ./bin
EXECUTABLE := $(BUILD_DIR)/main

SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC))

.PHONY: all clean cleanall

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(BUILD_DIR)/*.o

cleanall: clean
	rm -f $(EXECUTABLE)

default: all
