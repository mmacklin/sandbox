SHELL = /bin/sh
CC = g++
AR = ar

CFLAGS = -g -Wall -I.. -D_DEBUG 
ARFLAGS = cq 

TARGET  = core.a 
SOURCES = $(wildcard *.cpp) 
HEADERS = $(wildcard *.h)

OBJECTS = $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS) Makefile
	$(AR) $(ARFLAGS) $(TARGET) $(OBJECTS) 

clean:
	-rm -f $(OBJECTS)
	-rm -f $(TARGET)

%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<

run: $(TARGET)
	./$(TARGET)

.PHONY : all clean
