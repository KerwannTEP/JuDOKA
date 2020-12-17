# julia Cut.jl --parallel no --lmax 5
##################################################

##################################################

aMeasure = 10.0 # Semi-major axis (in mpc) along which the coefficients Djj are computed

##################################################
using HDF5 # To have access to .hf5 files
########################################
include("../sources/Main.jl") # Loading the main code
########################################
jminMeasure, jmaxMeasure = 0.01,0.999 # Range in j where the Djj are computed
nbjMeasure = 30 # Number of j for which the Djj are computed
tabjMeasure = exp.(range(log(jminMeasure),length=nbjMeasure,log(jmaxMeasure)))

const tabDRRjj = zeros(Float64,nbjMeasure) # Values of the DRR_jj coefficients on the (a,j)-grid
const tabDNRjj = zeros(Float64,nbjMeasure) # Values of the DRR_jj coefficients on the (a,j)-grid

########################################

function tabDRRjj!()
    if (PARALLEL == "yes") # Computation is made in parallel
        Threads.@threads for ij=1:nbjMeasure # Loop over the angular momenta
            jMeasure = tabjMeasure[ij] # Current angular momentum
            if jMeasure<jlc(aMeasure) # Do not compute DRR within the loss cone
                tabDRRjj[ij] = 0.0
            else
                tabSMARes_parallel = zeros(Float64,2,2) # Container of tabSMARes for the current thread 
                IntTable_K_parallel = IntTable_create!()
                
                DRR = DRR_jj(aMeasure,jMeasure,IntTable_K_parallel,tabSMARes_parallel)
                tabDRRjj[ij] = DRR # Computing DRR_jj with the thread's containers
            end
        end
    else # Computation is not made in parallel
        for ij=1:nbjMeasure # Loop over the angular momenta
            jMeasure = tabjMeasure[ij] # Current angular momentum
            if jMeasure<jlc(aMeasure) # Do not compute DRR within the loss cone
                tabDRRjj[ij] = 0.0
            else
                DRR = DRR_jj(aMeasure,jMeasure)
                tabDRRjj[ij] = DRR # Computing DRR_jj with the global containers
            end
        end
    end
end

function tabDNRjj!()
    if (PARALLEL == "yes") # Computation is made in parallel
        Threads.@threads for ij=1:nbjMeasure # Loop over the angular momenta
            jMeasure = tabjMeasure[ij] # Current angular momentum
            if jMeasure<jlc(aMeasure) # Do not compute DNR within the loss cone
                tabDNRjj[ij] = 0.0
            else
                DNR = DNR_jj(aMeasure,jMeasure)
                tabDNRjj[ij] = DNR # Computing DNR_jj with the thread's containers
            end
        end
    else # Computation is not made in parallel
        for ij=1:nbjMeasure # Loop over the angular momenta
            jMeasure = tabjMeasure[ij] # Current angular momentum
            if jMeasure<jlc(aMeasure) # Do not compute DNR within the loss cone
                tabDNRjj[ij] = 0.0
            else
                DNR = DNR_jj(aMeasure,jMeasure)
                tabDNRjj[ij] = DNR # Computing DNR_jj with the global containers
            end
        end
    end
end

########################################
namefile = "../data/Dump_Diffusion_Coefficients_Cut.hf5"
########################################
# Function that writes to .hf5 files
########################################
function writedump!(namefile)
    file = h5open(namefile,"w") # Opening the file
    write(file,"tabj",tabjMeasure) # Dumping the list of angular momenta
    write(file,"tabDRRjj",tabDRRjj) # Dumping the DRRjj(a,j) for the current used value of lmax
    write(file,"tabDNRjj",tabDNRjj) # Dumping the DNRjj(a,j) for the current used value of lmax
    write(file,"aMeasure",aMeasure) # Dumping aMeasure
    write(file,"jmin",jminMeasure) # Dumping jminMeasure
    write(file,"jmax",jmaxMeasure) # Dumping jmaxMeasure
    write(file,"nbj",nbjMeasure) # Dumping nbjMeasure
    write(file,"jlc",jlc(aMeasure)) # Dumping jlc(a)
    close(file) # Closing the file
end

########################################

@time tabDRRjj!()
@time tabDNRjj!()

########################################
writedump!(namefile) # Dumping the computed table
