# determine OS: Linux vs. Windows (mingw)
SYS := $(shell gcc -dumpmachine)

ifneq (, $(findstring mingw, $(SYS)))
	# Windows MinGW
	CC			= C:/Rtools/mingw64/bin/gcc.exe
	CFLAGS		= -std=gnu99 -O3
	R_C_HEADER	= C:/"Program Files"/R/R-4.0.2/include
	DLL_BLAS	= C:/"Program Files"/R/R-4.0.2/bin/x64/Rblas.dll
	DLL_LAPACK	= C:/"Program Files"/R/R-4.0.2/bin/x64/Rlapack.dll
	DLL_R		= C:/"Program Files"/R/R-4.0.2/bin/x64/R.dll
else
	# Linux
	CC			= gcc
	CFLAGS		= -g -std=gnu99 -O2
	R_C_HEADER	= /usr/share/R/include/
	DLL_BLAS	= /usr/lib/x86_64-linux-gnu/libblas.so
	DLL_LAPACK	= /usr/lib/x86_64-linux-gnu/liblapack.so
	DLL_R		= /usr/lib/R/lib/libR.so
endif

# link
ifneq (, $(findstring mingw, $(SYS)))
# Windows MinGW
robsurvey: robsurvey.o wquantile.o mallows.o zeroin.o huberm.o \
	trimmedwinsorized.o fitwls.o psifunctions.o robsurvey_error.o
	$(CC) -shared -o $@.dll $@.o $(DLL_BLAS) $(DLL_LAPACK) $(DLL_R) \
		wquantile.o mallows.o zeroin.o huberm.o \
		trimmedwinsorized.o fitwls.o psifunctions.o robsurvey_error.o
else
# Linux

endif

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
trimmedwinsorized.o: trimmedwinsorized.c
	$(CC) -Wall -pedantic -I $(R_C_HEADER) -c $*.c $(CFLAGS)

# compile
fitwls.o: fitwls.c
	$(CC) -Wall -pedantic -I $(R_C_HEADER) -c $*.c $(CFLAGS)

# compile
psifunctions.o: psifunctions.c
	$(CC) -Wall -pedantic -I $(R_C_HEADER) -c $*.c $(CFLAGS)

# compile with -fpic flag otherwise the array of strings is not relocatable
robsurvey_error.o: robsurvey_error.c
	$(CC) -fpic -Wall -pedantic -I $(R_C_HEADER) -c $*.c $(CFLAGS)

# housekeeping (using make clean)
.PHONY: clean
clean:
	rm robsurvey.o wquantile.o zeroin.o huberm.o mallows.o \
	trimmedwinsorized.o fitwls.o psifunctions.o robsurvey_error.o
