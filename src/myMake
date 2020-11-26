# some macros
CC		= C:/Rtools/mingw64/bin/gcc.exe
CFLAGS		= -std=gnu99 -O3
R_C_HEADER	= C:/"Program Files"/R/R-4.0.2/include
BLAS_LIBS	= C:/"Program Files"/R/R-4.0.2/bin/x64/Rblas.dll 
LAPACK_LIBS	= C:/"Program Files"/R/R-4.0.2/bin/x64/Rlapack.dll 
DLL_R		= C:/"Program Files"/R/R-4.0.2/bin/x64/R.dll

# link
robsurvey: robsurvey.o wquantile.o mallows.o zeroin.o huberm.o cbacon.o \
	trimmedwinsorized.o 
	$(CC) -shared -o $@.dll $@.o $(BLAS_LIBS) $(LAPACK_LIBS) $(DLL_R) \
		wquantile.o mallows.o zeroin.o huberm.o cbacon.o \
		trimmedwinsorized.o

# compile
robsurvey.o: robsurvey.c 
	$(CC) -Wall -pedantic -I $(R_C_HEADER) -c $*.c $(CFLAGS) 

# compile
wquantile.o: wquantile.c 
	$(CC) -Wall -pedantic -I $(R_C_HEADER) -c $*.c $(CFLAGS) 

# compile
mallows.o: mallows.c 
	$(CC) -Wall -pedantic -I $(R_C_HEADER) -c $*.c $(CFLAGS) 

# compile
zeroin.o: zeroin.c 
	$(CC) -Wall -pedantic -I $(R_C_HEADER) -c $*.c $(CFLAGS) 

# compile
huberm.o: huberm.c 
	$(CC) -Wall -pedantic -I $(R_C_HEADER) -c $*.c $(CFLAGS) 

# compile
cbacon.o: cbacon.c 
	$(CC) -Wall -pedantic -I $(R_C_HEADER) -c $*.c $(CFLAGS) 

# compile
trimmedwinsorized.o: trimmedwinsorized.c 
	$(CC) -Wall -pedantic -I $(R_C_HEADER) -c $*.c $(CFLAGS) 


# housekeeping (using make clean)
.PHONY: clean
clean: 
	rm robsurvey.o wquantile.o zeroin.o huberm.o mallows.o cbacon.o \
	trimmedwinsorized.o 
