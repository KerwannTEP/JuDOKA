# JuDOKA

Julia Diffusion of Orbits for Keplerian Actions

## INSTALLATION

Install Julia by following the instruction at `https://julialang.org/downloads/platform/`.

## PACKAGES

1) Open the terminal and type

```
$ julia
```

2) Install the following packages:

- `StaticArrays`
- `HypergeometricFunctions`
- `BenchmarkTools`
- `HDF5`
- `DelimitedFiles`
- `Interpolations`
- `ArgParse`
- `SpecialFunctions`
- `Distributions`
- `Plots`

For example, to install the package StaticArrays, one needs to run the command 
```
julia> import Pkg; Pkg.add("StaticArrays") 
```

To exit the Julia terminal, type the command
```
julia> exit()
```
### !! WARNING !!

**DO NOT INTERRUPT THE DOWNLOADING OF THE PACKAGES !!!!**

## MODIFY THE BACKGROUND FAMILIES

To modify the background of the S-Cluster, open the folder `code/background`
and open the file `inputBath.txt`.
	
Do not change the line 1, which serves as legend.
Starting from line 2, enter the values for each family, one family per line.
Enter the four values, separated by a tabulation.

The numbers describe:
- `gamma`: Cusp power of the family.
- `mass_star[Msun]`: Mass of the family's individuals, which can be stars or BHs (in solar masses).
- `a0[mpc]`: scale factor of the family (in mpc).
- `Mtot(<a0)`: Enclosed mass of the family within the radius a0 (scale factor given previously).

## PLOT A DIFFUSION COEFFICIENT a-cut MAP IN ORBITAL SPACE

To compute a mapping of the diffusion coefficients DRR and DNR at fixed 
semi-major axis, and depending on reduced angular momentum, one needs to open 
`code/tests/Cut.jl`.

Then, one needs to write for which sma `a` one wants to compute DRRjj(a,j). 
Once this is done, save the file and run the command 

```
$ julia Cut.jl
```

within the folder `code/tests`. If one wants to run this with parallelization,
one needs to run the following commands (supposing one is using bash)

```
$ export JULIA_NUM_THREADS=4
$ export JULIA_CPU_THREADS=4
$ julia -p 4 Cut.jl --parallel yes --lmax 10
```

where 4 is the number of parallelized threads. One can check the number of 
threads by opening the Julia terminal and by running the command

```
julia> Threads.nthreads()
```

The resulting file will be created in the folder `code/data` under the name 
`Dump_Diffusion_Coefficients_Cut.hf5`.

Once one has recovered the `Dump_Diffusion_Coefficients_Cut.hf5` file using the latter
method, one can plot the a-cut in orbital space of the relevant coefficients.

### PLOT USING JULIA

Go to the folder `code/tests/Julia` and run the file `PlotCut.jl` using the command

```
$ julia PlotCut.jl
```
    
A log-log plot of the diffusion coefficients which have just been computed is
displayed on a window. To terminate the programm, close the window and type
the **`Enter`** key on the terminal.

Once this is done, the figure is recovered as a PNG file of the name `DjjCut.png` 
in the folder `code/graphs/Julia`.

### PLOT USING MATHEMATICA

Go to the folder `code/tests/Scripts` and evaluate the Mathematica script
`Diffusion_Coefficients_Cut.m` using the command

```
$ math -script Diffusion_Coefficients_Cut.m 
```

Once done, we recover the plot as a PNG file of the name `DjjCut.png` 
in the folder "code/graphs/Mathematica".


## PLOT A DIFFUSION COEFFICIENT MAP IN ORBITAL SPACE

To compute an orbital space map of the SRR diffusion coefficients DRRjj and DNRjj,
one needs to access the `code/tests` folder and run the following command in 
the terminal:

```
$ julia Mapping.jl
```

If one wants to run this with parallelization, one needs to run the following 
commands (supposing one is using bash)

```
$ export JULIA_NUM_THREADS=4
$ export JULIA_CPU_THREADS=4
$ julia -p 4 Mapping.jl --parallel yes --lmax 10
```
	
where 4 is the number of parallelized threads. One can check the number of 
threads by opening the Julia terminal and by running the command

```
julia> Threads.nthreads()
```

The resulting file will be created in the folder `code/data` under the name 
`Dump_Diffusion_Coefficients.hf5`.

### PLOT USING JULIA

Go to the folder `code/tests/Julia` and run the file `PlotContours.jl` using the command

```
$ julia PlotContours.jl
```
    
A log-log plot of the diffusion coefficients which have just been computed is
displayed on a window. To terminate the programm, close the window and type
the **`Enter`** key on the terminal.

Once this is done, the figure is recovered as a PNG file of the name `Djj.png` 
in the folder `code/graphs/Julia`.

### PLOT USING MATHEMATICA

Go to the folder `code/tests/Scripts` and evaluate the Mathematica script
`Diffusion_Coefficients.m` using the command

```
$ math -script Diffusion_Coefficients.m 
```

Once done, we recover the plot as a PNG file of the name `Djj.png` 
in the folder `code/graphs/Mathematica`.

