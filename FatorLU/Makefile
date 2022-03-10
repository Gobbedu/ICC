    CC     = gcc -g
    CFLAGS = -g -std=c11
    LFLAGS = -Wl,-rpath -Wl, -lm

    PROG = testaSL
    OBJS = utils.o SistLinear.o EliminacaoGauss.o Refinamento.o

.PHONY: clean purge all

%.o: %.c
	$(CC) -c $(CFLAGS) $<

$(PROG):  $(OBJS) testaSL.o
	$(CC) -o $@ $^ $(LFLAGS)

clean:
	@rm -f *~ *.bak

purge:   clean
	@rm -rf bin obj *.layout *.depend
	@rm -f *.o core a.out
	@rm -f $(PROG)
