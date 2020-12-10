#using DelimitedFiles # To be able to load .txt files

##################################################
# Function that returns the input data from the cluster
# The array has the shape {nbStar,4},
# and is stored in the order {nameStar,a(arcsec),e,age}
# Used to get nbStar and to fill in
# taba,tabe,tabAge
##################################################

const d_BH = 8.32*10^6 # Distance to Sgr A* (in mpc)
const conversion_factor = d_BH * 2*pi/(360*60*60)

function arcsec_to_mpc(a::Float64)
    return a * conversion_factor
end

function ecc_to_j(ecc::Float64)
    return sqrt(1.0-ecc^2)
end



function uncertainty_ecc_to_j(ecc::Float64,de::Float64)
    j = ecc_to_j(ecc)
    dj = ecc/j * de
    return dj
end

function get_dataCluster()
    dataCluster, headerCluster = readdlm(INPUTCLUSTER,header=true) # Reading the data file with the header
    return dataCluster # Output the entire data
end

const dataCluster = get_dataCluster() # Reading the input data from the Cluster
const nbCluster = size(dataCluster)[1] # Number of stars in cluster

#################################################
# Getting the relevant physical quantities of the S-cluster
##################################################

const tabaCluster = SVector{nbCluster,Float64}([arcsec_to_mpc(dataCluster[iCluster,2]) for iCluster=1:nbCluster]) # The sma a is stored in column 2 (mpc)
const tabjCluster = SVector{nbCluster,Float64}([ecc_to_j(dataCluster[iCluster,3]) for iCluster=1:nbCluster]) # The eccentricity e is stored in column 3
const tabTCluster = SVector{nbCluster,Float64}([dataCluster[iCluster,4] for iCluster=1:nbCluster]) # The main-sequence age is stored in column 4 (in Myr)

const tabdjCluster = SVector{nbCluster,Float64}([uncertainty_ecc_to_j(dataCluster[iCluster,3],dataCluster[iCluster,5]) for iCluster=1:nbCluster]) # The eccentricity uncertainty de is stored in column 5
const tabdTpCluster = SVector{nbCluster,Float64}([dataCluster[iCluster,6] for iCluster=1:nbCluster]) # The age up-uncertainty dT_p is stored in column 6
const tabdTmCluster = SVector{nbCluster,Float64}([dataCluster[iCluster,7] for iCluster=1:nbCluster]) # The age down-uncertainty dT_m is stored in column 7
