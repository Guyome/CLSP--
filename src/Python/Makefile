CC=gcc
CFLAGS=-Wall -ansi -pedantic -O3 -fPIC
LDFLAGS=-shared -lblitz -lboost_python -lpython2.6 -llapack -lipopt -lpthread
SRCDIR=./
TESTDIR=$(SRCDIR)/test
LIB=heurclsp.so
PYTHONDIR=/usr/include/python2.6
COINDIR=/usr/include/coin
SRC= $(wildcard $(SRCDIR)/*.cpp)
OBJ= $(SRC:.cpp=.o)

all: $(LIB)

$(LIB): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)
	cp $(SRCDIR)/$(LIB) $(TESTDIR)

%.o: %.cpp
	$(CC) -o $@ -c $< $(CFLAGS) -I $(PYTHONDIR) -I $(COINDIR)

clean:
	rm -rf $(SRCDIR)/*.o *.so
	rm -rf $(TESTDIR)/*.so


