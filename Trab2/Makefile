
    CC     = gcc -g
    CFLAGS = -g -std=c11 -O3 -DLIKWID_PERFMON -I/home/soft/likwid
    LFLAGS = -Wl,-rpath -Wl, -lm -lmatheval

    PROG = newtonPC
    # OBJS = NewtonPadrao.o NewtonModificado.o NewtonInexato.o SistNlinear.o SistLinear.o utils.o
    OBJS = NewtonPadrao.o Rosenbrock.o NewtonInexato.o SistNlinear.o SistLinear.o utils.o

.PHONY: clean purge all

%.o: %.c
	$(CC) -c $(CFLAGS) $<

$(PROG):  $(OBJS) newtonPC.o
	$(CC) -o $@ $^ $(LFLAGS)

clear: 
	@rm $(OBJS)

purge:
	@rm -rf bin obj *.layout *.depend
	@rm -f *.o core a.out
	@rm -f $(PROG)
