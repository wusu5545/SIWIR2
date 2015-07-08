CXX = g++
CXXFLAGS = -Wall -pedantic -std=c++11

INC_DIR = include
SRC_DIR = src
OBJ_DIR = object
OUT_DIR = vtk
 
TARGET = mdsim
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
HEDS = $(wildcard $(INC_DIR)/*.h)
OBJS = $(SRCS:$(SRC_DIR)%.cpp=$(OBJ_DIR)%.o)

argu  = input/mini.par  input/mini.dat
argu1 = input/params_50x50.dat
argu2 = input/params_bonus.dat

REF_DIR = referenceOutputs/vtk_50x50
REF_SRC =$(wildcard $(REF_DIR)/*.vtk)

VTK_DIR = vtk


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

cleanvtk:
	@$(RM) -rf $(OUT_DIR)/*.vtk

clean:
	@$(RM) -rf deps.dat $(OBJ_DIR)/*.o $(TARGET) $(OBJ_DIR) $(OUT_DIR)/*.vtk
	
rebuild: clean $(TARGET)

compare:
	@for i in $(notdir $(REF_SRC)); do \
	echo "Compare $$i"; \
	cmp -s $(VTK_DIR)/$$i $(REF_DIR)/$$i; \
	RETVAL=$$?; \
	if [ $$RETVAL -eq 0 ]; then \
		echo "SAME"; \
	else \
		echo "NOT SAME"; \
	fi \
	done;

run: $(TARGET) cleanvtk
	./$(TARGET) $(argu)
	
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

