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
## WARNING !!

**DO NOT INTERRUPT THE DOWNLOADING OF THE PACKAGES !!!!**

## MODIFY THE BACKGROUND FAMILIES

To modify the background of the S-Cluster, open the folder "code/background"
and open the file "inputBath.txt".
	
Do not change the line 1, which serves as legend.
Starting from line 2, enter the values for each family, one family per line.
Enter the four values, separated by a tabulation.

The numbers describe:
- gamma: Cusp power of the family.
- mass_star[Msun]: Mass of the family's individuals, which can be stars or BHs (in solar masses).
- a0[mpc]: scale factor of the family (in mpc).
- Mtot: Enclosed mass of the family within the radius a0 (scale factor given previously).

	
