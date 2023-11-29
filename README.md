# 1DDiffusion

C and MATLAB code to compute the numerical solution of the 1D Diffusion equation using Crank-Nicolson differencing on a periodic domain.

## Build & Run C code

1. Change directory into the /C folder and run the Makefile


```sh
cd ./C
make
```

running `make` creates an executable called `diff`.

2. Run the executable for default parameters (these can be manually changed in the diffusion.c file)


```sh
./diff
``` 

this will run the code and save the solution data in `Data/C-Data` directory

## Plot Data

You can analyse the solution data by running the plot.py script

```sh
python3 plot.py
``` 

this script computes the $l_2$ norm of the numerical solution and the exact solution over time.
