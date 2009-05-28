CC=gcc
CFLAGS=-W -Wall -ansi -pedantic -O3
LDFLAGS=-shared
LIBS=-lbliz
LIB=libclsp.so
SRC= $(wildcard src/*.cpp)
OBJ= $(SRC:.cpp=.o)

all: $(LIB)

$(LIB): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CC) -o $@ -c $< $(CFLAGS) $(LIBS)

clean:
	rm -rf src/*.o
