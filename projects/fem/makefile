SHELL = /bin/sh
CC = g++

CFLAGS = -g -O3 -ffast-math -Wall -I"../.." -D_DEBUG
LDFLAGS = -g -framework GLUT -framework OpenGL -framework Cocoa

TARGET  = fem 
SOURCES = $(wildcard *.cpp) ../../core/platform.cpp ../../core/shader.cpp ../../core/maths.cpp ../../core/tga.cpp 
HEADERS = $(wildcard *.h)

OBJECTS = $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $(TARGET) 

clean:
	-rm -f $(OBJECTS)
	-rm -f $(TARGET)

%.o: %.cpp $(HEADERS) makefile
	$(CC) $(CFLAGS) -c -o $@ $<

run: $(TARGET)
	./$(TARGET)

.PHONY : all clean
