# Makefile
#! /bin/sh

DEBUG=no

CC = gcc
ifeq ($(DEBUG),yes)
	CFLAGS= -Wextra -g
else
	CFLAGS= -O3
endif

#macos
LDFLAGS = -I/usr/local/include/graphviz/
#linux
#LDFLAGS += -L/usr/lib/graphviz/
#LDFLAGS += -I/usr/include/graphviz/

LDFLAGS += -lgvc -lcgraph -lcdt -lm


SRC = $(wildcard *.c)
OBJ= $(SRC:.cpp=.o)
EXEC = neattest

all: $(EXEC)
ifeq ($(DEBUG),yes)
	@echo "Génération en mode debug"
else
	@echo "Génération en mode release"
endif

$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@ 

.PHONY: clean mrproper

clean:
	rm -f *.o

mrproper: clean
	rm -f $(EXEC)

