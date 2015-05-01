SHELL := /bin/bash
NAME = figure


# set executable name
EXECNAME    = main

# default compiler settings
CC          =  clang++
OPT         = -O3
IFLAGS     = -I/usr/lib/gcc/x86_64-pc-cygwin/4.9.2/include -I/usr/lib/gcc/x86_64-pc-cygwin/4.9.2/include/c++ -I/usr/lib/gcc/x86_64-pc-cygwin/4.9.2/include/c++/x86_64-pc-cygwin/ -std=c++11 -stdlib=libc++
LDFLAGS     = 

SRC         = *.cc
OBJS        = $*(SRC).o

# compilation for runs
all:
	$(CC) $(IFLAGS) -g $(OPT) $(SRC) -o $(EXECNAME) $(LDFLAGS) -v -Wno-error=unknown-warning
	/bin/rm -rf *.o

# clean up
clean:
	rm -rf *.o $(EXECNAME)
