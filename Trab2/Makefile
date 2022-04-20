
    CC     = gcc -g -O3 -mavx -march=native
	LIKWIDF= -DLIKWID_PERFMON -I${LIKWID_INCLUDE} -L${LIKWID_LIB} -llikwid
    LFLAGS = -Wl,-rpath -Wl, -lm -lmatheval

    PROG = newtonPC
    # OBJS = NewtonPadrao.o NewtonModificado.o NewtonInexato.o SistNlinear.o SistLinear.o utils.o
    OBJS = NewtonPadrao.o Rosenbrock.o NewtonInexato.o SistNlinear.o SistLinear.o utils.o

.PHONY: clean purge all

# gcc -DLIKWID_PERFMON -I${LIKWID_INCLUDE} -c <prog.c>
# gcc -o <prog> <prog.o> -L${LIKWID_LIB} -llikwid

%.o: %.c
	$(CC) $(LIKWIDF) -c $<

$(PROG):  $(OBJS) newtonPC.o
	$(CC) $(LIKWIDF) -o $@ $^ $(LFLAGS)

clear: 
	@rm $(OBJS)

purge:
	@rm -rf bin obj *.layout *.depend
	@rm -f *.o core a.out
	@rm -f $(PROG)
