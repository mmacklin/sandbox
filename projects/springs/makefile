SHELL = /bin/sh
CC = g++

CFLAGS = -g -m32 -Wall -I"../.." -D_DEBUG
LDFLAGS = -g -m32 -framework GLUT -framework OpenGL -framework Cocoa

TARGET  = Springs
SOURCES = $(wildcard *.cpp) ../../core/platform.cpp ../../core/shader.cpp
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
