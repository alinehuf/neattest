# Makefile
#! /bin/sh

DEBUG=yes
SYSTEM=mac

CC = gcc
ifeq ($(DEBUG),yes)
	CFLAGS= -Wextra -g
else
	CFLAGS= -O3
endif

ifeq ($(SYSTEM),win)
	LIBS = -lGL -lGLU -lglut -lm -ljpeg       # win
	CFLAGS += -DWIN
else
	ifeq ($(SYSTEM),linux)
		LIBS = -lGL -lGLU -lglut -lm -ljpeg       # linux
		CFLAGS += -DLINUX
	else
		LIBS = -framework Glut -framework OpenGL  # mac
    CFLAGS += -DMAC
	endif
endif
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
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@ 

.PHONY: clean mrproper

clean:
	rm -f *.o

rmproper: clean
	rm -f $(EXEC)

