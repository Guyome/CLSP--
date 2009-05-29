CC=gcc
CFLAGS=-W -Wall -ansi -pedantic
LDFLAGS=-lblitz
EXEC=example
SRC=example.cpp ../HeurClsp.cpp
all: $(EXEC)

$(EXEC): $(SRC)
	$(CC) -o $@ $^ $(LDFLAGS) $(CFLAGS)

clean:
	rm -rf *.o
