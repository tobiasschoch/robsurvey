# some macros
CC		= C:/Rtools/mingw_64/bin/gcc.exe
CFLAGS		= -std=c11 -O0
R_C_HEADER	= C:/"Program Files"/R/R-3.5.1/include
BLAS_LIBS	= C:/"Program Files"/R/R-3.5.1/bin/x64/Rblas.dll 
LAPACK_LIBS	= C:/"Program Files"/R/R-3.5.1/bin/x64/Rlapack.dll 
DLL_R		= C:/"Program Files"/R/R-3.5.1/bin/x64/R.dll

# link
robsurvey: robsurvey.o wquantile.o mallows.o zeroin.o 
	$(CC) -shared -o $@.dll $@.o $(BLAS_LIBS) $(LAPACK_LIBS) $(DLL_R) \
 wquantile.o mallows.o zeroin.o

# compile
robsurvey.o: robsurvey.c 
	$(CC) -Wall -I $(R_C_HEADER) -c $*.c $(CFLAGS) 

# compile
wquantile.o: wquantile.c 
	$(CC) -Wall -I $(R_C_HEADER) -c $*.c $(CFLAGS) 

# compile
mallows.o: mallows.c 
	$(CC) -Wall -I $(R_C_HEADER) -c $*.c $(CFLAGS) 

# compile
zeroin.o: zeroin.c 
	$(CC) -Wall -I $(R_C_HEADER) -c $*.c $(CFLAGS) 

# housekeeping (using make clean)
.PHONY: clean
clean: 
	rm robsurvey.o wquantile.o zeroin.o 
