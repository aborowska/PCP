mex -v -largeArrayDims -ldl CFLAGS="\$CFLAGS -std=c99" myfun.c -lmwblas -lmwlapack


mex -v -ldl CFLAGS="\$CFLAGS -std=c99" C_posterior_garch11_mex.c
mex -v -ldl CFLAGS="\$CFLAGS -std=c99" C_addparam_posterior_garch11_mex.c
mex -v -ldl CFLAGS="\$CFLAGS -std=c99" posterior_garch11_mex.c
