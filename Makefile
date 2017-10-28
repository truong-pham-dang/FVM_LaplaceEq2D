CF         = gfortran
FFLAGS     = -O3 -ffixed-line-length-160 -fbounds-check -g -fdollar-ok
LD         = gfortran
LDFLAGS    = 
PREPROC    = 

OBJS =  main.o \
        gauss_2.o \
        matinv.o  \


.SUFFIXES: .o .f90 .f
.f90.o:
	$(LD) -c $(FFLAGS) $<
.f.o:
	$(LD) -c $(FFLAGS) $<

Lap :$(OBJS) 
	$(LD) $(LDFLAGS) -o $@ $(OBJS)

clean :
	rm -f Lap *.o core *.mod

