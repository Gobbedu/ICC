    CC     = gcc -g -O3 -mavx -march=native
    # LFLAGS = -Wl,-rpath -Wl, -lm -lmatheval 
	LIKWID = -L${LIKWID_LIB} -llikwid

    PROG = newtonPC
    OBJS = NewtonPadrao.o NewtonInexato.o Rosenbrock.o SistNlinear.o SistLinear.o utils.o
    CODE = NewtonPadrao.c NewtonInexato.c Rosenbrock.c SistNlinear.c SistLinear.c utils.c

.PHONY: clean purge all

%.o: %.c
	$(CC) -DLIKWID_PERFMON -I${LIKWID_INCLUDE} -c $<

$(PROG):  $(OBJS) newtonPC.o
	$(CC) -o $@ $^ $(LIKWID) $(LFLAGS)

# para compilar sem likwid
local:
	$(CC) -c $(CODE) newtonPC.c
	$(CC) -o $(PROG) newtonPC.o $(OBJS) $(LFLAGS)


clear: 
	@rm $(OBJS)

purge:
	@rm -rf bin obj *.layout *.depend
	@rm -f *.o core a.out
	@rm -f $(PROG)