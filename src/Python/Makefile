CC=gcc
CFLAGS=-Wall -W -O3 -fPIC -pipe
LDFLAGS=-shared -lblitz -lboost_python -lpython2.6 -llapack -lipopt -lpthread
SRCDIR=./
TESTDIR=$(SRCDIR)/test
LIB=heurclsp.so
INCLUDE=/usr/include/python2.6
SRC= $(wildcard $(SRCDIR)/*.cpp)
OBJ= $(SRC:.cpp=.o)

all: $(LIB)

$(LIB): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)
	cp $(SRCDIR)/$(LIB) $(TESTDIR)

%.o: %.cpp
	$(CC) -o $@ -c $< $(CFLAGS) -I $(INCLUDE)

clean:
	rm -rf $(SRCDIR)/*.o *.so
	rm -rf $(TESTDIR)/*.so


