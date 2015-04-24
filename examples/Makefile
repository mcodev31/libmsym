LELEMENTS = msym msym_error context symmetry symmetrize linalg point_group character_table geometry orbital\
		permutation symop equivalence_set elements

SRC = example.c
OBJS = example.o
LSRC = $(foreach I,$(LELEMENTS),../src/$I.c)
LOBJS = $(foreach I,$(LELEMENTS),../src/$I.o)



CC = gcc
EXE = msym_example
CFLAGS = -g -O2 -flto -std=c99
INC = ../src
LIBS = -lm

all:    $(OBJS) $(LOBJS)
	$(CC) -o $(EXE) $(CFLAGS) $(OBJS) $(LOBJS) $(LIBS)

%.o: %.c
	$(CC) ${CFLAGS} -I $(INC) -o $@ -c $<

clean:  
	$(RM) -f $(OBJS) $(LOBJS)

