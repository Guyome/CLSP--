CC=gcc
CFLAGS=-W -Wall -ansi -pedantic -O3 -fPIC
LDFLAGS=-shared -lblitz -lboost_python -lpython2.6
SRCDIR=src
LIB=libclsp.so
PYLIB=heurclsp.so
PYDIR=/usr/include/python2.6
TEST=example
TESTDIR=$(SRCDIR)/test
SRC= $(wildcard $(SRCDIR)/*.cpp)
OBJ= $(SRC:.cpp=.o)

all: $(LIB) $(TEST) $(PYLIB)

$(LIB): $(SRCDIR)/HeurClsp.o
	$(CC) -o $@ $^ $(LDFLAGS)

$(PYLIB): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

$(TEST):
	@(cd $(TESTDIR) && $(MAKE) $@)

%.o: %.cpp
	$(CC) -o $@ -c $< $(CFLAGS) -I $(PYDIR)

clean:
	rm -rf src/*.o *.so
	@(cd $(TESTDIR) && $(MAKE) $@)

