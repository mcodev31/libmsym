LELEMENTS = msym msym_error context symmetry symmetrize linalg point_group character_table geometry orbital\
		permutation symop equivalence_set elements vibration subspace

SRC = example.c msym_example.c tex_character_table.c
OBJDIR = ../obj
BINDIR = ../bin
EXAMPLE_OBJS = $(OBJDIR)/example.o $(OBJDIR)/msym_example.o
TEX_OBJS = $(OBJDIR)/tex_character_table.o 
LSRC = $(foreach I,$(LELEMENTS),../src/$I.c)
LOBJS = $(foreach I,$(LELEMENTS),$(OBJDIR)/$I.o)



CC = gcc
EXAMPLE_EXE = $(BINDIR)/msym_example
TEX_EXE = $(BINDIR)/tex_character_table
CFLAGS = -g -O2 -flto -std=c99 -pedantic -Wall
INC = ../src
LIBS = -lm

all:    $(OBJDIR) $(BINDIR) $(EXAMPLE_OBJS) $(TEX_OBJS) $(LOBJS)
	$(CC) -o $(EXAMPLE_EXE) $(CFLAGS) $(EXAMPLE_OBJS) $(LOBJS) $(LIBS)
	$(CC) -o $(TEX_EXE) $(CFLAGS) $(TEX_OBJS) $(LOBJS) $(LIBS)

$(OBJDIR) :
	mkdir -p $(OBJDIR)

$(BINDIR) :
	mkdir -p $(BINDIR)

../obj/example.o: example.c
	$(CC) ${CFLAGS} -I $(INC) -o $@ -c $<

../obj/msym_example.o: msym_example.c
	$(CC) ${CFLAGS} -I $(INC) -o $@ -c $<

../obj/tex_character_table.o: tex_character_table.c
	$(CC) ${CFLAGS} -I $(INC) -o $@ -c $<

../obj/%.o: ../src/%.c
	$(CC) ${CFLAGS} -I $(INC) -o $@ -c $<

clean:  
	$(RM) -rf $(OBJDIR) $(BINDIR)

