FORTRAN  =  pgf90
OPTS     = 
LOADER   =  pgf90
LOADOPTS =  -O4 -fast -mp -mcmodel=medium -o exe
LIBRARY  = -lblas 

FILES  =  main.o zcopy.o zdscal.o zscal.o zaxpy.o zdotc.o makeS.o makeRHS.o AS.o BS.o CS.o DS.o ZS.o KS.o MS.o NS.o PS.o U_Mult.o bicgstabU.o B_Mult.o S_solve.o params.o



OBJS   = 

exe : Makefile $(FILES) 
	$(LOADER) $(LOADOPTS) $(FILES) $(OBJS) $(LIBRARY) 
.f.o:
	pgf90 -O4 -fast -mp -mcmodel=medium -c -g $*.f
clean:
	rm -f *.o *.inp core a.out exe
 
