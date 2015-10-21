CXX = g++
CXXFLAGS = -Wall -std=c++11 -pedantic

INC_DIR = include
SRC_DIR = src
OBJ_DIR = object
 
TARGET = waveguide
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
HEDS = $(wildcard $(INC_DIR)/*.h)
OBJS = $(SRCS:$(SRC_DIR)%.cpp=$(OBJ_DIR)%.o)

argu = 0.01 10e-10 2

.PHONY: mkdir clean all deps push

all: $(TARGET)

deps: $(SRCS) Makefile
	@echo -n "building dependencies..."
	@rm -rf deps.dat
	@for i in $(SRCS); do \
		echo "$(OBJ_DIR)/" | tr -d '\n' >> deps.dat; \
		gcc -std=c++11 -MM $$i -I $(INC_DIR) >> deps.dat; \
		echo "" >> deps.dat; \
	done;
	@echo " Success!"

-include deps.dat

mkdir:
ifeq ($(wildcard $(OBJ_DIR)),)
	@mkdir $(OBJ_DIR)
endif

$(TARGET): mkdir deps $(OBJS)
	@echo -n "compiling $@..."
	@$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)
	@echo " Success!"

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp Makefile
	@echo -n "compiling $@..."
	@$(CXX) $(CXXFLAGS) -I $(INC_DIR) -o $@ -c $<
	@echo " Success!"

clean:
	@$(RM) -rf deps.dat $(OBJ_DIR)/*.o $(TARGET) $(OBJ_DIR)
	
rebuild: clean $(TARGET)

run: $(TARGET)
	./$(TARGET) $(argu)
	gnuplot gnuplotScript


MODULE = default

new:
ifeq ($(MODULE),default)
	@echo 'usage: make new MODULE=<newmodule>'
else
ifeq ($(MODULE),main)
	@echo '\n\nint main(int argc, char* argv[]){\n\treturn 0;\n}' >> $(SRC_DIR)/$(MODULE).cpp
else
	@echo '#include "$(MODULE).h"' > $(SRC_DIR)/$(MODULE).cpp
	@echo -n '#ifndef ' > $(INC_DIR)/$(MODULE).h
	@echo -n '$(MODULE)' | tr '[:lower:]' '[:upper:]' >> $(INC_DIR)/$(MODULE).h
	@echo '_H' >> $(INC_DIR)/$(MODULE).h
	@echo -n '#define ' >> $(INC_DIR)/$(MODULE).h
	@echo -n '$(MODULE)' | tr '[:lower:]' '[:upper:]' >> $(INC_DIR)/$(MODULE).h
	@echo '_H' >> $(INC_DIR)/$(MODULE).h
	@echo '' >> $(INC_DIR)/$(MODULE).h
	@echo '' >> $(INC_DIR)/$(MODULE).h
	@echo '' >> $(INC_DIR)/$(MODULE).h
	@echo -n '#endif//' >> $(INC_DIR)/$(MODULE).h
	@echo -n '$(MODULE)' | tr '[:lower:]' '[:upper:]' >> $(INC_DIR)/$(MODULE).h
	@echo '_H' >> $(INC_DIR)/$(MODULE).h
endif
	@echo 'created module "$(MODULE)"!'
endif

MESSAGE = default

push:
ifeq ($(MESSAGE),default)
	@echo 'usage: make push MESSAGE=<message>'
else
	@echo 'committing...'
	@git commit -a -m "$(MESSAGE)"
	@echo 'successfully committed!'
	@echo 'pushing...'
	@git push
	@echo 'successfully pushed!'
endif

