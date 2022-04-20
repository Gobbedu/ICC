# # gcc -DLIKWID_PERFMON -I${LIKWID_INCLUDE} -c <prog.c>
# # gcc -o <prog> <prog.o> -L${LIKWID_LIB} -llikwid

    CC     = gcc -g -O3 -mavx -march=native
    LFLAGS = -Wl,-rpath -Wl, -lm -lmatheval -L${LIKWID_LIB} -llikwid

    PROG = newtonPC
    OBJS = NewtonPadrao.o NewtonInexato.o Rosenbrock.o SistNlinear.o SistLinear.o utils.o

.PHONY: clean purge all

%.o: %.c
	$(CC) -DLIKWID_PERFMON -I${LIKWID_INCLUDE} -c $<

$(PROG):  $(OBJS) newtonPC.o
	$(CC) -o $@ $^ $(LFLAGS)

clear: 
	@rm $(OBJS)

purge:
	@rm -rf bin obj *.layout *.depend
	@rm -f *.o core a.out
	@rm -f $(PROG)