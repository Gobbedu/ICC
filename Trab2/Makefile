
    CC     = gcc -g -O3 -mavx -march=native
    CFLAGS = -g -std=c11 
	LIKWIDF= -DLIKWID_PERFMON -I /home/soft/likwid/include -L /home/soft/likwid/lib -llikwid
    LFLAGS = -Wl,-rpath -Wl, -lm -lmatheval

    PROG = newtonPC
    # OBJS = NewtonPadrao.o NewtonModificado.o NewtonInexato.o SistNlinear.o SistLinear.o utils.o
    OBJS = NewtonPadrao.o Rosenbrock.o NewtonInexato.o SistNlinear.o SistLinear.o utils.o

.PHONY: clean purge all

%.o: %.c
	$(CC) $(LIKWIDF) -c $(CFLAGS)  $<

$(PROG):  $(OBJS) newtonPC.o
	$(CC) $(LIKWIDF) -o $@ $^ $(LFLAGS)

clear: 
	@rm $(OBJS)

purge:
	@rm -rf bin obj *.layout *.depend
	@rm -f *.o core a.out
	@rm -f $(PROG)
