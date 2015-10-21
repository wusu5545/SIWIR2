VPATH = src: include: test

CC = g++ 
CFLAGS =  -c -O3 -Wall -Winline -Wshadow -ansi -I./myColsamm/Source

INCLUDES = 

SOURCES = waveguide.cpp FiniteElement.cpp
DEPS = FiniteElement.h
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = waveguide

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ 

.cpp.o:
	$(CC) $(CFLAGS)  $< -o $@ 
	
clean:
	-rm -f $(OBJECTS)
