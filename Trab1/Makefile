    CC     = gcc -g
    CFLAGS = -g -std=c11
    LFLAGS = -Wl,-rpath -Wl, -lm -lmatheval

    PROG = testaSNL
    OBJS = utils.o SistLinear.o EliminacaoGauss.o FatoracaoLU.o SistNlinear.o

.PHONY: clean purge all

%.o: %.c
	$(CC) -c $(CFLAGS) $<

$(PROG):  $(OBJS) testaSNL.o
	$(CC) -o $@ $^ $(LFLAGS)

clear: 
	@rm $(OBJS)

clean:
	@rm -f *~ *.bak

purge:   clean
	@rm -rf bin obj *.layout *.depend
	@rm -f *.o core a.out
	@rm -f $(PROG)
