##################################################
# Computes a random walk
# Command: julia RandomWalk.jl --inputBath ../../../code/background/inputBath.txt --inputCluster ../../../code/background/inputCluster.txt

# Recover the graph in the folder code/graphs
##################################################

using Random, Distributions # To be able to use the normal distribution
using DelimitedFiles # To be able to load .txt files
using Plots # To be able to plot data
using LaTeXStrings
using HDF5

include("../sources/Main.jl") # Loading the main code

#################################################
# Parameters used in the computation of the diffusion coefficients
##################################################

const aWalk = 10.0 # Semi-major axis where the random walk happens

const nbjMeasure = 100 # Number of j for which the coefficients are computed
const jminMeasure, jmaxMeasure = 0.001,1.0 # Range in j where the coefficients are computed
const tabjMeasure = range(jminMeasure,length=nbjMeasure,jmaxMeasure)
const arrayDj_Djj_serial = zeros(Float64,2,nbjMeasure) # Values of the D_j and D_jj coefficients on the (a,j)-grid
const arrayTjj_serial = zeros(Float64,nbjMeasure) # Values of the D_j and D_jj coefficients on the (a,j)-grid

##################################################
# Parameters of the stochastic walk
# Numbers of stars, of runs, initial conditions IC, etc...
##################################################

const IC = 0.6 # Initial condition is a Dirac of "nb_stars" stars at j=IC, for t=0.0 Myr here
const timeEnd = 10.0 # End time (in Myr)

#################################################
# Initializing the coefficients' tables
##################################################

function init_arrayDj_Djj!(arrayDj_Djj=arrayDj_Djj_serial)
    for ij=1:nbjMeasure
        arrayDj_Djj[1,ij] = 0.0
        arrayDj_Djj[2,ij] = 0.0
        arrayTjj_serial[ij] = 0.0
    end
end

#################################################
# Functions to compute the diffusion coefficients
##################################################

const jlca = jlc(aWalk)

const tabDjj_nr_rr = zeros(Float64, 2, nbjMeasure)

function arrayDj_Djj!(arrayDj_Djj_arg=arrayDj_Djj_serial,arrayTjj_arg=arrayTjj_serial,Table_dKdJ=IntTable_dKdJ_serial,Table_dKdJp=IntTable_dKdJp_serial,tabSMARes=tabSMARes_serial)
    init_arrayDj_Djj!(arrayDj_Djj_arg) # Making sure that the table is initially set to 0
    ijmin = 1
    for ij=1:nbjMeasure # Loop over the elements of the j
        j = tabjMeasure[ij]
        if j<jlca # Do not compute DRR within the loss cone
            arrayTjj_serial[ij] = 0.0
            ijmin = ij
        elseif j<1.0
            DRR = DRR_j_jj(aWalk,j,Table_dKdJ,Table_dKdJp,tabSMARes) # Computing DRR
            DNR = DNR_aj(aWalk,j) # Computing DNR
            arrayDj_Djj_arg[1,ij] = DRR[1] + DNR[3]
            arrayDj_Djj_arg[2,ij] = DRR[2] + DNR[4]

            # Timestep criterion is taken from Shapiro & Marchant (1978), applied for both NR and SRR, following Bar-Or & Alexander (2016)
            tnr = min(0.01/DNR[4], 0.16*(1.005-j)^2/DNR[4])
            trr = min(0.01/DRR[2], 0.16*(1.005-j)^2/DRR[2])
            arrayTjj_serial[ij] = min(tnr, trr)

            tabDjj_nr_rr[1, ij] = DNR[4]
            tabDjj_nr_rr[2, ij] = DRR[2]
        else
            arrayDj_Djj_arg[1,ij] = 0.0
            arrayDj_Djj_arg[2,ij] = 0.0
            arrayTjj_serial[ij] = arrayTjj_serial[ij-1]
        end
    end  
    return ijmin+1
end

println("Computing the diffusion coefficients...")

@time ijmin = arrayDj_Djj!()

p = plot(tabjMeasure[ijmin:nbjMeasure-1], [arrayDj_Djj_serial[2,ijmin:nbjMeasure-1] tabDjj_nr_rr[1,ijmin:nbjMeasure-1] tabDjj_nr_rr[2,ijmin:nbjMeasure-1]],
        yaxis=:log10, xaxis=:log10, 
        xticks=10.0 .^ (-3:1:0), xminorticks=10,
        yticks=10.0 .^ (-6:1:2), yminorticks=10,
        frame=:box, label=["Total" "NR" "RR"],
        xlabel=L"j", ylabel=L"D_{jj}",
        title="Diffusion coefficients jj",
        legend=:topleft)

display(p)
readline()

#println( arrayDj_Djj[1,:])
#println( arrayDj_Djj[2,:])
#################################################
# Functions to Interpolate the diffusion coefficients
##################################################

function get_IntDj(arrayDj_Djj=arrayDj_Djj_serial) # Returns the interpolation function of Dj at fixed a_{iCluster}
    intDj = Interpolations.scale(interpolate(arrayDj_Djj[1,:], BSpline(Linear())),tabjMeasure) # Constructing the interpolation function for j-> Dj(j,a_Cluster)
    return intDj
end

function get_IntDjj(arrayDj_Djj=arrayDj_Djj_serial) # Returns the interpolation function of Djj at fixed a_{iCluster}
    intDjj = Interpolations.scale(interpolate(arrayDj_Djj[2,:], BSpline(Linear())),tabjMeasure) # Constructing the interpolation function for j-> Djj(j,a_Cluster)
    return intDjj
end

function get_IntTjj(arrayTjj=arrayTjj_serial) # Returns the interpolation function of Djj at fixed a_{iCluster}
    intTjj = Interpolations.scale(interpolate(arrayTjj, BSpline(Linear())),tabjMeasure) # Constructing the interpolation function for j-> Djj(j,a_Cluster)
    return intTjj
end


##################################################
# Driving the stochastic process
##################################################

function step_j(j::Float64,alea::Float64,stepT::Float64,intDj,intDjj)
    Dj = intDj(j)
    Djj = intDjj(j)
    return Dj*stepT + sqrt(Djj*stepT) * alea 
end

function evolve(t::Float64,j::Float64,intDj,intDjj) # Drives the stochastic process until T_star
    alea = rand(Normal())
    jc = jlca # Angular momentum at loss-con 
    tau = intTjj(j)
    println("(t, dt, j) = ", (t, tau, j)) 
    if j<=jc
        return (timeEnd,j)
    else
        stepJ = step_j(j,alea,tau,intDj,intDjj)
        if ((j + stepJ<1.0) && (j + stepJ>0.0))
            return (t+tau, j+stepJ)
        elseif (j + stepJ>0.0)
            return (t+tau, 2.0-(j+stepJ)) # Elastic rebound of the star-particle against bound j=1          
       else 
            return (t+tau, abs(j+stepJ)) # Elastic rebound of the star-particle against bound j=0
        end
    end
end    

println("Interpolating the diffusion coefficients...")

intDj = get_IntDj()
intDjj = get_IntDjj()
intTjj = get_IntTjj()

tStar = 0.0
jStar = IC

tab = [tStar jStar]


namefile = "RandomWalk.txt"

file = open(namefile, "a")
writedlm(file, tab)

println("Running the random walk...")

while (tStar < timeEnd)
    tstart = tStar
    global tStar, jStar = evolve(tStar, jStar, intDj, intDjj)
    tabs = [tStar jStar]
    writedlm(file, tabs)   

    dt = tStar-tstart
    if (dt == 0.0)
        break
    end 

    if (jStar <= jlca)
        break 
    end
    
    
end

close(file)

data = readdlm(namefile,header=false)
nb = size(data)[1]

rm(namefile)

const tabt = zeros(Float64,nb)
const tabj = zeros(Float64,nb)
const tabdt = zeros(Float64,nb-1)

for it=1:nb
    tabt[it] = data[it,1]
    tabj[it] = data[it,2]
end

for it=2:nb 
    tabdt[it-1] = data[it,1]-data[it-1,1]
end

log10dtmin = floor(Int64, log10(minimum(tabdt)))
log10dtmax = floor(Int64, log10(maximum(tabdt))) + 1


p = plot(tabt, tabj, legend=false,
        xlabel="t [Myr]", ylabel=L"j",
        title="Eccentricity relaxation",
        yticks=0:0.2:1, yminorticks=2,
        ylims=(0, 1),
        xticks=0:2:10, xminorticks=2,
        xlims=(0, timeEnd),
        frame=:box)

plot!(p, [0, timeEnd], [jlca, jlca], label=:false, 
            linewidth=1,
            linestyle=:dash,
            alpha=1,
            color=:red)

annotate!(p, 0.9*timeEnd, jlca+0.02, text("Loss cone", :red, 7,  rotation = 0))

savefig(p,"../graphs/Julia/RandomWalk_j.png") # Saves the figure
display(p) # Display plot
readline() # Plot window stays open until we press "Enter"



p = plot(tabt[2:end], tabdt, legend=false,
        xlabel="t [Myr]", ylabel="dt [Myr]",
        title="Timestep",
        yticks=10.0 .^ (-5:1:0), yminorticks=10,
        ylims=(10.0^(log10dtmin), 10.0^(log10dtmax)),
        yaxis=:log10,
        xticks=0:2:10, xminorticks=2,
        xlims=(0, timeEnd),
        frame=:box)

savefig(p,"../graphs/Julia/RandomWalk_dt.png") # Saves the figure
display(p) # Display plot
readline() # Plot window stays open until we press "Enter"


newfile = "../data/RandomWalk.hf5"

function writedump!(newfile)
    file = h5open(newfile,"w")
    write(file,"tabt",tabt)
    write(file,"tabdt",tabdt)
    write(file,"tabj",tabj)
end

writedump!(newfile)
