VPATH = src: include: test

CC = g++ 
CFLAGS =  -c -O3 -Wall -Winline -Wshadow -ansi -fopenmp
LDFLAGS = -lgomp -lpthread

INCLUDES = -I/

SOURCES = main.cpp mgm.cpp functions.cpp solveMG.cpp fmg.cpp
DEPS = FUNCTIONS.h
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = mgsolve

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS)  $< -o $@
	
clean:
	-rm -f $(OBJECTS)
