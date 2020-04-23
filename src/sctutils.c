void print_double(int n, int p, double *x);
void print_integer(int n, int p, int *x);


/*****************************************************************************\
|*  print matrix							     *|
|*                                                                           *|
|*  PARAMETERS                                                               *|
|*    n		  dimension                                                  *|
|*    p		  dimension                                                  *|
|*    x		  array[n, p]	                                             *|
|*                                                                           *|
|*  NOTE: for vectors, take n = 1                                            *|
|*                                                                           *|
\*****************************************************************************/
void print_double(int n, int p, double *x){
   for (int i = 0; i < n; i++ ) {
      for(int j = 0; j < p; j++ ) Rprintf("%10.3f\t", x[i + j * n]);
   Rprintf("\n");
   }
   Rprintf("\n");
}

void print_integer(int n, int p, int *x){
   for (int i = 0; i < n; i++ ) {
      for(int j = 0; j < p; j++ ) Rprintf("%d\t", x[i + j * n]);
   Rprintf("\n");
   }
   Rprintf("\n");
}




