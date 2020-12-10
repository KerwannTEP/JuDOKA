# julia Mapping.jl --parallel no --lmax 5
##################################################
# Parameters of the grid in (a,j)--space
# This is the grid where the DRR & DNR coefficients are computed
# Evalue a between aminMeasure (mpc) and amaxMeasure (mpc)
##################################################
using HDF5 # To have access to .hf5 files
########################################
include("../sources/Main.jl") # Loading the main code
########################################
jminMeasure, jmaxMeasure = 0.001,1.0 # Range in j where the Djj are computed
aminMeasure, amaxMeasure = 0.1,1000.0 # Range in j where the Djj are computed
nbjMeasure = 300 # Number of j for which the Djj are computed
nbaMeasure = 100 # Number of a for which the Djj are computed
nbajGrid = nbjMeasure*nbaMeasure # Number of (a,j) for which the Djj are computed
tabjMeasure = exp.(range(log(jminMeasure),length=nbjMeasure,log(jmaxMeasure)))
tabaMeasure = exp.(range(log(aminMeasure),length=nbaMeasure,log(amaxMeasure)))

const tabajGrid  = zeros(Float64,2,nbajGrid) # Location (a,j) of the grid points where the diffusion coefficients are computed
const tabDRRjjGrid = zeros(Float64,nbajGrid) # Values of the DRR_jj coefficients on the (a,j)-grid
const tabDNRjjGrid = zeros(Float64,nbajGrid) # Values of the DNR_jj coefficients on the (a,j)-grid
const tabDjjGrid = zeros(Float64,nbajGrid) # Values of the total D_jj coefficients on the (a,j)-grid
const tabDNRJJRedGrid = zeros(Float64,nbajGrid) # Values of the DNR_jj coefficient$

########################################

function tabajGrid!()
    index = 1
    for ia=1:nbaMeasure
    aMeasure = tabaMeasure[ia]
        for ij=1:nbjMeasure
            jMeasure = tabjMeasure[ij]
            tabajGrid[1,index], tabajGrid[2,index] = aMeasure, jMeasure
            index += 1
        end
    end
end

function init_tabDRRGrid!()
    for iGrid=1:nbajGrid
        tabDRRjjGrid[iGrid] = 0.0
    end
end

function init_tabDNRGrid!()
    for iGrid=1:nbajGrid
        tabDNRjjGrid[iGrid] = 0.0
        
        tabDNRJJRedGrid[iGrid] = 0.0
    end
end

function init_tabDGrid!()
    for iGrid=1:nbajGrid
        tabDjjGrid[iGrid] = 0.0
    end
end

########################################

function tabDRRGrid!()
    init_tabDRRGrid!() # Making sure that the grid is initially set to 0
    if (PARALLEL == "yes") # Computation is made in parallel
        Threads.@threads for iGrid=1:nbajGrid # Loop over the elements of the (a,j)-grid
            a, j = tabajGrid[1,iGrid], tabajGrid[2,iGrid] # Current (a,j) location
            if j<=jlc(a) # Do not compute DRR within the loss cone
                tabDRRjjGrid[iGrid] = 0.0
            elseif j<1.0
                tabSMARes_parallel = zeros(Float64,4) # Container of tabSMARes for the current thread 
                Table_K_parallel = Table_create!()
                
                DRR = DRR_jj(a,j,Table_K_parallel,tabSMARes_parallel)
                tabDRRjjGrid[iGrid] = DRR # Computing DRR_jj with the thread's containers
            else
                tabDRRjjGrid[iGrid] = 0.0
            end
        end
    else # Computation is not made in parallel
        for iGrid=1:nbajGrid # Loop over the elements of the (a,j)-grid
            a, j = tabajGrid[1,iGrid], tabajGrid[2,iGrid] # Current (a,j) location
            if j<=jlc(a) # Do not compute DRR within the loss cone
                tabDRRjjGrid[iGrid] = 0.0
            elseif j<1.0
                DRR = DRR_jj(a,j)
                tabDRRjjGrid[iGrid] = DRR # Computing DRR_jj with the global containers
            else
                tabDRRjjGrid[iGrid] = 0.0
            end
        end
    end
end

function tabDNRGrid!()
    init_tabDNRGrid!() # Making sure that the grid is initially set to 0
    if (PARALLEL == "yes") # Computation is made in parallel
        Threads.@threads for iGrid=1:nbajGrid # Loop over the elements of the (a,j)-grid
            a, j = tabajGrid[1,iGrid], tabajGrid[2,iGrid] # Current (a,j) location
#            if j<=jlc(a) # Do not compute DNR within the loss cone
 #               tabDNRjjGrid[iGrid] = 0.0
                
  #              tabDNRJJRedGrid[iGrid] = 0.0
   #         elseif j<1.0   
            if j<1.0 
                DNR = DNR_jj(a,j)
                tabDNRjjGrid[iGrid] = DNR # Computing DNR_jj with the thread's containers

                tabDNRJJRedGrid[iGrid] = DNR_JJ_red(a,j)
            else
                tabDNRjjGrid[iGrid] = 0.0
                tabDNRJJRedGrid[iGrid] = 0.0
            end
        end
    else # Computation is not made in parallel
        for iGrid=1:nbajGrid # Loop over the elements of the (a,j)-grid
            a, j = tabajGrid[1,iGrid], tabajGrid[2,iGrid] # Current (a,j) location
 #            if j<=jlc(a) # Do not compute DNR within the loss cone
 #               tabDNRjjGrid[iGrid] = 0.0
                
  #              tabDNRJJRedGrid[iGrid] = 0.0
   #         elseif j<1.0   
            if j<1.0 
                DNR = DNR_jj(a,j)
                tabDNRjjGrid[iGrid] = DNR # Computing DNR_jj with the global containers
            
                tabDNRJJRedGrid[iGrid] = DNR_JJ_red(a,j)
            else
                tabDNRjjGrid[iGrid] = 0.0
                tabDNRJJRedGrid[iGrid] = 0.0
            end
        end
    end
end

function tabDGrid!()
    init_tabDGrid!() # Making sure that the grid is initially set to 0
    if (PARALLEL == "yes") # Computation is made in parallel
        Threads.@threads for iGrid=1:nbajGrid # Loop over the elements of the (a,j)-grid
            a, j = tabajGrid[1,iGrid], tabajGrid[2,iGrid] # Current (a,j) location
  #          if j<=jlc(a) # Do not compute D_jj within the loss cone
   #             tabDjjGrid[iGrid] = 0.0
    #        elseif j<1.0 
            if j<1.0    
                DNR = tabDNRjjGrid[iGrid]
                DRR = tabDRRjjGrid[iGrid]
                tabDjjGrid[iGrid] = DNR + DRR # Computing D_jj with the thread's containers
            else
                tabDjjGrid[iGrid] = 0.0
            end
        end
    else # Computation is not made in parallel
        for iGrid=1:nbajGrid # Loop over the elements of the (a,j)-grid
            a, j = tabajGrid[1,iGrid], tabajGrid[2,iGrid] # Current (a,j) location
        #    if j<=jlc(a) # Do not compute D_jj within the loss cone
       #         tabDjjGrid[iGrid] = 0.0
       #     elseif j<1.0  
            if j<1.0   
                DNR = tabDNRjjGrid[iGrid]
                DRR = tabDRRjjGrid[iGrid]
                tabDjjGrid[iGrid] = DNR + DRR # Computing D_jj with the thread's containers
            else
                tabDjjGrid[iGrid] = 0.0
            end
        end
    end
end

########################################
namefile = "../data/Dump_Diffusion_Coefficients.hf5"
########################################
# Function that writes to .hf5 files
########################################
function writedump!(namefile)
    file = h5open(namefile,"w") # Opening the file
    write(file,"tabaj",tabajGrid) # Dumping the grid of (a,j)
    write(file,"tabj",tabjMeasure) # Dumping the grid of j
    write(file,"taba",tabaMeasure) # Dumping the grid of a
    write(file,"tabDRRjj",tabDRRjjGrid) # Dumping the DRRjj(a,j) for the current used value of lmax
    write(file,"tabDNRjj",tabDNRjjGrid) # Dumping the DNRjj(a,j) for the current used value of lmax
    write(file,"tabDjj",tabDjjGrid) # Dumping the total Djj(a,j) for the current used value of lmax
    write(file,"amin",aminMeasure) # Dumping aminMeasure
    write(file,"amax",amaxMeasure) # Dumping amaxMeasure
    write(file,"jmin",jminMeasure) # Dumping jminMeasure
    write(file,"jmax",jmaxMeasure) # Dumping jmaxMeasure
    write(file,"nbj",nbjMeasure) # Dumping nbjMeasure
    write(file,"nba",nbaMeasure) # Dumping nbaMeasure
    
    
    write(file,"tabDNRJJRed",tabDNRJJRedGrid)
    
    write(file,"G",G)
    write(file,"cvel",cvel)
    write(file,"mBH",mBH)
    
    
    
    close(file) # Closing the file
end

########################################

tabajGrid!()

@time tabDRRGrid!()
@time tabDNRGrid!()
@time tabDGrid!()

########################################
writedump!(namefile) # Dumping the computed table
