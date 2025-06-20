# Compiler and flags
CC = gcc
NVCC = nvcc
CFLAGS = -O3
NVCCFLAGS = -O3 -arch=sm_60  # Adjust the architecture flag based on your GPU
LIBFLAGS = -lm -lcudart      # Include the CUDA runtime library

OBJDIR = obj

# Object files
_OBJ = args.o data.o setup.o vtk.o lid-cavity.o
OBJ = $(patsubst %,$(OBJDIR)/%,$(_OBJ))

.PHONY: directories clean all

# Default target
all: directories lid-cavity

# Rule for object files
$(OBJDIR)/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

$(OBJDIR)/%.o: %.cu
	$(NVCC) -c -o $@ $< $(NVCCFLAGS)

# Link the executable
lid-cavity: $(OBJ)
	$(NVCC) -o $@ $^ $(CFLAGS) $(LIBFLAGS)

# Create directories
directories: $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

# Clean up
clean:
	rm -Rf $(OBJDIR)
	rm -f lid-cavity
