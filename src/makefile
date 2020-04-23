IDIR =../hdr
CC=gcc
CFLAGS=-I$(IDIR) -ggdb

ODIR=obj
LDIR =../lib

LIBS=-lm -lgsl -lgslcblas 

_DEPS = fem2d.h gauss_legendre.h config.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = curing_2d_nonlinear.o nodes_connectivity.o material_properties.o local_elements.o vector_products.o gauss_seidel.o global_matrix.o curing_cycle.o matrix_properties.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

_EXTOBJ = gauss_legendre.o 
EXTOBJ = $(patsubst %,$(ODIR)/%,$(_EXTOBJ))


$(ODIR)/%.o: %.c
		$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: $(LDIR)/%.c
		$(CC) -c -o $@ $< $(CFLAGS)

curing_2d_nonlinear: $(OBJ) $(EXTOBJ)
		gcc -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
		rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
