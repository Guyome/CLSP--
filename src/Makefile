LIB=libclsp.so
PYLIB=heurclsp.so
TEST=example
LIBDIR=C++
PYLIBDIR=Python
TESTDIR=$(LIBDIR)/test

all: $(LIB) $(TEST) $(PYLIB)

$(LIB): 
	@(cd $(LIBDIR) && $(MAKE) $@)
	mv $(LIBDIR)/$(LIB) $(LIB)

$(PYLIB):
	@(cd $(PYLIBDIR) && $(MAKE) $@)
	mv $(PYLIBDIR)/$(PYLIB) $(PYLIB)

$(TEST):
	@(cd $(TESTDIR) && $(MAKE) $@)

clean:
	rm -rf $(SRCDIR)/*.o *.so
	@(cd $(LIBDIR) && $(MAKE) $@)
	@(cd $(PYLIBDIR) && $(MAKE) $@)

