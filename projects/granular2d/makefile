SHELL = /bin/sh
CC = g++

CFLAGS = -g -m32 -Wall -I"../.." -O3 -ffast-math -msse -msse2 -msse3 -mfpmath=sse 
LDFLAGS = -g -m32 -framework GLUT -framework OpenGL -framework Cocoa 

TARGET  = granular2d

SOURCES = $(wildcard *.cpp) ../../core/platform.cpp ../../core/shader.cpp ../../core/maths.cpp ../../core/mesh.cpp ../../core/aabbtree.cpp ../../core/tga.cpp
HEADERS = $(wildcard *.h)

OBJECTS = $(SOURCES:.cpp=.o) 

all: $(TARGET)

$(TARGET): $(OBJECTS) makefile
	$(CC) $(LDFLAGS) $(OBJECTS) -o $(TARGET) 

clean:
	-rm -f $(OBJECTS)
	-rm -f $(TARGET)

%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<


run: $(TARGET)
	./$(TARGET)

.PHONY : all clean
