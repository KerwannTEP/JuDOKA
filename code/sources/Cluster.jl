#using DelimitedFiles # To be able to load .txt files

##################################################
##################################################

"""
    conversion_factor_arcsec_to_sma

Multiplicative conversion factor from arcseconds to mpc.
"""
const conversion_factor_arcsec_to_sma = d_BH * 2*pi/(360*60*60)

"""
    arcsec_to_mpc(a::Float64)

Converts a semi-major axis given in arcseconds to an semi-major axis given in mpc.
"""
function arcsec_to_mpc(a::Float64)
    return a * conversion_factor_arcsec_to_sma
end

"""
    get_dataCluster()

Recovers the S-cluster data from INPUTCLUSTER.
"""
function get_dataCluster()
    dataCluster, headerCluster = readdlm(INPUTCLUSTER,header=true) # Reading the data file with the header
    return dataCluster # Output the entire data
end

"""
    dataCluster

Table of S-cluster data
"""
const dataCluster = get_dataCluster() # Reading the input data from the Cluster

"""
    nbCluster

Number of S-cluster stars used in the code
"""
const nbCluster = size(dataCluster)[1] # Number of stars in cluster

#################################################
# Getting the relevant physical quantities of the S-cluster
##################################################

"""
    tabaCluster

Table of S-cluster SMAs
"""
const tabaCluster = SVector{nbCluster,Float64}([arcsec_to_mpc(dataCluster[iCluster,2]) for iCluster=1:nbCluster]) # The sma a is stored in column 2 (mpc)

"""
    tabjCluster

Table of S-cluster reduced angular momenta.
"""
const tabjCluster = SVector{nbCluster,Float64}([ecc_to_j(dataCluster[iCluster,3]) for iCluster=1:nbCluster]) # The eccentricity e is stored in column 3

"""
    tabTCluster

Table of S-cluster main-sequence ages.
"""
const tabTCluster = SVector{nbCluster,Float64}([dataCluster[iCluster,4] for iCluster=1:nbCluster]) # The main-sequence age is stored in column 4 (in Myr)

"""
    tabdjCluster

Table of S-cluster SMA uncertainties.
"""
const tabdjCluster = SVector{nbCluster,Float64}([uncertainty_ecc_to_j(dataCluster[iCluster,3],dataCluster[iCluster,5]) for iCluster=1:nbCluster]) # The eccentricity uncertainty de is stored in column 5

"""
    tabdTpCluster

Table of S-cluster upper main-sequence age uncertainties.
"""
const tabdTpCluster = SVector{nbCluster,Float64}([dataCluster[iCluster,6] for iCluster=1:nbCluster]) # The age up-uncertainty dT_p is stored in column 6

"""
    tabdTmCluster

Table of S-cluster lower main-sequence age uncertainties.
"""
const tabdTmCluster = SVector{nbCluster,Float64}([dataCluster[iCluster,7] for iCluster=1:nbCluster]) # The age down-uncertainty dT_m is stored in column 7