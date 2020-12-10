##################################################
# Use this file to compute the diffusion coefficients
# Use the cluster model given in the "sources" folder

# Use the following command to get truncated quantities up to lmax = 10
# julia Compute.jl --lmax 10
##################################################

include("../sources/Main.jl")

a = 10.0 # Semi-major axis (in mpc)
j = 0.6 # Reduced angular momentum J/Jc (between 0 and 1)

println("----------------------------------------------")
println("Computing both DRRj and DRRjj at the same time")
@time djBoth, djjBoth = DRR_j_jj(a,j)
println("DRRj(a,j)  = ",djBoth)
println("DRRjj(a,j) = ",djjBoth)

println("----------------------------------------------")

println("Computing DRRjj only")
@time djj_only = DRR_jj(a,j)
println("DRRjj(a,j) = ",djj_only)
println("----------------------------------------------")
