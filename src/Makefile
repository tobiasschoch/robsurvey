CC         = C:/Rtools/mingw64/bin/gcc.exe
CFLAGS     = -std=gnu99 -O3

R_DIR      = C:/"Program Files"/R/R-4.0.2
R_C_HEADER = $(R_DIR)/include
DLL_BLAS   = $(R_DIR)/bin/x64/Rblas.dll
DLL_LAPACK = $(R_DIR)/bin/x64/Rlapack.dll
DLL_R      = $(R_DIR)/bin/x64/R.dll

# Link
OBJ = wquantile.o mallows.o zeroin.o huberm.o trimmedwinsorized.o \
		psifunctions.o robsurvey_error.o robsurvey.o

robsurvey.dll : $(OBJ)
	$(CC) -shared -o $@ $(DLL_BLAS) $(DLL_LAPACK) $(DLL_R) $(OBJ)

# Compile
%.o : %.c
	$(CC) -Wall -pedantic -I $(R_C_HEADER) -c -o $@ $< $(CFLAGS)

# Housekeeping
.PHONY: clean
clean:
	rm $(OBJ)
