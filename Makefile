CC := g++
SRCDIR := src
BUILDDIR := build
TESTDIR := test
GTESTDIR := /home/mackendy/apps/gtest-1.7.0
TARGET := bin/graph

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
TEST_SRCS := $(shell find $(TESTDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
TEST_OBJS := $(patsubst $(TESTDIR)/%,$(BUILDDIR)/%,$(TEST_SRCS:.$(SRCEXT)=.o))
CPPFLAGS = -std=c++14 -O2 -Wall -g -DDEBUG
LDFLAGS =
INC := -Iinclude/ -I/usr/include/boost
LIB := -lboost_program_options -lpthread
TESTINC := $(INC) -I$(GTESTDIR)/include
TESTLIB := $(LIB) -L$(GTESTDIR)/lib/.libs -lgtest -lgtest_main

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	$(CC) $^ -o $(TARGET) $(LIB)

$(TARGET)_test: $(TEST_OBJS) $(filter-out $(BUILDDIR)/main.o, $(OBJECTS))
	$(CC) $^ -o $(TARGET)_test $(TESTLIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	$(CC) $(CPPFLAGS) $(INC) -c -o $@ $<

$(BUILDDIR)/%.o: $(TESTDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	$(CC) $(CPPFLAGS) $(TESTINC) -c -o $@ $<

test: $(TARGET)_test

clean:
	@echo " Cleaning...";
	$(RM) -r $(BUILDDIR) $(TARGET)

.PHONY: clean test
