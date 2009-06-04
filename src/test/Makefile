CC=gcc
CFLAGS=-W -Wall -ansi -pedantic
LDFLAGS=-lblitz -lboost_python -lpython2.6
INCDIR=/usr/include/python2.6
EXEC=example
SRC=example.cpp ../HeurClsp.cpp
all: $(EXEC)

$(EXEC): $(SRC)
	$(CC) -o $@ $^ $(LDFLAGS) $(CFLAGS) -I $(INCDIR)

clean:
	rm -rf *.o
