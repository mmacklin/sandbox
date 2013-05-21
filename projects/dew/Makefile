SHELL = /bin/sh
CC = g++

CFLAGS = -g -Wall -I"../.." 
LDFLAGS = -g -framework GLUT -framework OpenGL -framework Cocoa

TARGET  = Dew 
SOURCES = $(wildcard *.cpp) ../../core/platform.cpp ../../core/shader.cpp ../../Core/maths.cpp ../../Core/tga.cpp ../../Core/mesh.cpp dialog.mm
HEADERS = $(wildcard *.h)

OBJECTS = $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS) Makefile
	$(CC) $(LDFLAGS) $(OBJECTS) -o $(TARGET) 

clean:
	-rm -f $(OBJECTS)
	-rm -f $(TARGET)

%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<

run: $(TARGET)
	./$(TARGET)

.PHONY : all clean
