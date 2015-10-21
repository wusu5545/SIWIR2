VPATH = src: include: test

CC = g++ 
CFLAGS =  -c -O3 -Wall -ansi -pedantic

INCLUDES = 

SOURCES = lbm.cpp functions.cpp 
DEPS = FUNCTIONS.h GRID.h
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = lbm

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ 

.cpp.o:
	$(CC) $(CFLAGS)  $< -o $@ 
	
clean:
	-rm -f $(OBJECTS)
