SHELL := /bin/bash
NAME = figure


# set executable name
EXECNAME    = main

# default compiler settings
CC          =  g++
OPT         = -O3
LDFLAGS     = -lm

SRC         = *.cc
OBJS        = $*(SRC).o

# compilation for runs
all:
	$(CC) -g $(OPT) $(SRC) -o $(EXECNAME) $(LDFLAGS)
	/bin/rm -rf *.o

# clean up
clean:
	rm -rf *.o $(EXECNAME)
