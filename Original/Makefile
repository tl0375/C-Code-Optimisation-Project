CC=gcc
CFLAGS=-O3
LIBFLAGS=-lm

OBJDIR = obj

_OBJ = args.o data.o setup.o vtk.o lid-cavity.o
OBJ = $(patsubst %,$(OBJDIR)/%,$(_OBJ))

.PHONY: directories

all: directories lid-cavity

obj/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) 

lid-cavity: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBFLAGS) 

clean:
	rm -Rf $(OBJDIR)
	rm -f lid-cavity

directories: $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

