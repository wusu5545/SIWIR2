CXX = g++
CXXFLAGS = -std=c++11 -Wshadow -Wall -march=native -Winline -O3 -fopenmp #-ftree-vectorizer-verbose=2 #-g #-Wextra -Wshadow -Werror -unroll-loops -mavx -DNDEBUG 



INC_DIR = inc
SRC_DIR = src
OBJ_DIR = obj
 
TARGET = mgsolve
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
HEDS = $(wildcard $(INC_DIR)/*.h)
OBJS = $(SRCS:$(SRC_DIR)%.cpp=$(OBJ_DIR)%.o)

#Number of levels
LEVELS = 5
#Number of V-cycles
VNUMB = 10

BONISON = 0

.PHONY: mkdir clean all deps push vimConfig vimConfigRevert


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
	@echo " Success!\n"

obj/%.o: $(SRC_DIR)/%.cpp Makefile
	@echo -n "compiling $@..."
	@$(CXX) $(CXXFLAGS) -I $(INC_DIR) -o $@ -c $<
	@echo " Success!"

hsizeError: $(TARGET)
	./errorPlot.sh

clean:
	@$(RM) -rf deps.dat obj/*.o $(TARGET) init.dat solution.dat obj

rebuild: clean $(TARGET)

L = 11

run: $(TARGET)
	./$(TARGET) $(L)
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
	@echo '\n\n\n#endif' >> $(INC_DIR)/$(MODULE).h
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

vimConfig:
	@mv ~/.vimrc ~/.vimrc_old
	@cp .vimrc ~/.vimrc
	@echo "vim config set! Revert with make gitConfigRevert"

vimConfigRevert:
	@mv ~/.vimrc_old ~/.vimrc
