VPATH = src: include: test

CC = g++ 
CFLAGS =  -c -O3 -Wall -ansi -pedantic -std=c++11

INCLUDES = 

SOURCES = mdsim.cpp ParameterReader.cpp Setup.cpp
DEPS = ParameterReader.hpp Setup.hpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = mdsim

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ 

.cpp.o:
	$(CC) $(CFLAGS)  $< -o $@ 
	
clean:
	-rm -f $(OBJECTS)
