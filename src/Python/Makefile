CC=gcc
CFLAGS=-W -Wall -ansi -pedantic -O3 -fPIC
LDFLAGS=-shared -lblitz -lboost_python -lpython2.6
SRCDIR=./
TESTDIR=$(SRCDIR)/test
LIB=heurclsp.so
DIR=/usr/include/python2.6
SRC= $(wildcard $(SRCDIR)/*.cpp)
OBJ= $(SRC:.cpp=.o)

all: $(LIB)

$(LIB): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)
	cp $(SRCDIR)/$(LIB) $(TESTDIR)

%.o: %.cpp
	$(CC) -o $@ -c $< $(CFLAGS) -I $(DIR)

clean:
	rm -rf $(SRCDIR)/*.o *.so
	rm -rf $(TESTDIR)/*.so


