# makefile for ocean model.

.SUFFIXES: .o .h .f .F

FC = ifort -convert big_endian

.F.o:
	ifort -I. -c -O0 $<

OBJS = diag_month.o

sbc:	$(OBJS)
	$(FC) -O0 ${OBJS} -o sbc

$(OBJS):

clean:
	-rm -f *.o
