CC=gcc
CFLAGS=-W -Wall -ansi -pedantic -O3 -fPIC
LDFLAGS=-shared -lblitz
SRCDIR=src
LIB=libclsp.so
TEST=example
TESTDIR=$(SRCDIR)/test
SRC= $(wildcard $(SRCDIR)/*.cpp)
OBJ= $(SRC:.cpp=.o)

all: $(LIB) $(TEST)

$(LIB): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

$(TEST):
	@(cd $(TESTDIR) && $(MAKE) $@)

%.o: %.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf src/*.o
	@(cd $(TESTDIR) && $(MAKE) $@)

