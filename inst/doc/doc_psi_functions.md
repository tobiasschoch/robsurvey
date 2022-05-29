# Adding Support for Other $\psi$-Functions

Tobias Schoch

## 1 Introduction

The `robsurvey` package implements the Huber and Tukey (biweight)
$\psi$-functions. The functions are implemented in the C language,
see `src/psifunctions.c`. For the Huber $\psi$-function, the
*standard* function and an asymmetric $\psi$-function are implemented.

The functions are referenced by an integer value (in the C and R source code):

* `psi = 0`: Huber;
* `psi = 1`: asymmetric Huber;
* `psi = 2`: Tukey biweight.

For each type of $\psi$-function, the following three functions (in
the C language) must be defined:

* `psi`-function, $\psi(x)$, the actual $\psi$-function;
* `weight`-function, $w(x),$ associated with the $\psi$-function;
* `psi-prime`-function, the first derivative of the
  $\psi$-function, $\psi'(x)$.

The $\psi$-, $w$-, and $\psi'$-functions
have the same signature, which is shown here for a dummy function `foo()`.

```c
double foo(double x, const double k)
{
    # the code goes here
}
```

Argument `x` is the function argument and argument `k`
is the robustness tuning constant.

> **Limitations.** In this note, we consider only adding support
> for $\psi$-functions whose signature comply with the above dummy function.
> If you want to add functions that do not comply, you have to modify the
> existing code.

The method dispatch takes place in the functions (see
`src/psifunctions.c`):

* `get_wgt_function()`
* `get_psi_function()`
* `get_psi_prime_function()`

and is implemented with function pointers.

## 2 Adding another function

In order to add support for additional $\psi$-functions (which comply with
the above signature), follow these steps:

* add the C code of the new `psi`-, `weight`, and `psi-prime`-functions to the source;
* add an entry in the `switch` statement of the functions `get_wgt_function()`,
  `get_psi_function()`, and `get_psi_prime\_function()`;
* add/ modify the R-functions.
