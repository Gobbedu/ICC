    CC     = gcc -g

    CFLAGS = -I/usr/include
    LFLAGS = -lm -L/usr/lib -lmatheval

      PROG = labZero
      OBJS = utils.o \
             ZeroFuncao.o

.PHONY: limpa faxina clean purge distclean all

%.o: %.c %.h utils.h
	$(CC) -c $(CFLAGS) $<

$(PROG):  $(OBJS) $(PROG).o
	$(CC) -o $@ $^ $(LFLAGS)

limpa:
	@rm -f *~ *.bak

purge:   limpa
	@rm -f *.o core a.out
	@rm -f $(PROG)
