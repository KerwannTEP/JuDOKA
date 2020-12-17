include("Packages.jl") #import the packages
##################################################
include("Args.jl") # Parsing the command-line
##################################################
# General parameters
##################################################
"""
    PARALLEL
    
Determining if the code is run in parallel.
"""
const PARALLEL = parsed_args["parallel"]
if ((PARALLEL != "yes") && (PARALLEL != "no"))
    error("ERROR: UNKNOWN PARALLEL") # Unknown parallel procedure
end


"""
    lmax
    
Maximum harmonic number used in the multipole expansion of DRRjj.
"""
const lmax = parsed_args["lmax"]

"""
    nbK
    
Number of sampling points for geometric factor Klnnp.
"""
const nbK = parsed_args["nbK"]

"""
    nbResPoints
    
Number of points on the resonance lines.
"""
const nbResPoints = parsed_args["nbResPoints"]

"""
    mBH
    
Mass of the central BH, in solar masses (MSun).
"""
const mBH = parsed_args["mBH"]

"""
    rh
    
Influence radius of the BH, in milliparsecs (mpc).
We do not compute any resonance lines beyond that point.
"""
const rh = parsed_args["rh"]

"""
    INPUTBATH
    
Name of the input file for the bath.
"""
const INPUTBATH = parsed_args["inputBath"] 

"""
    INPUTCLUSTER

Name of the input file for the S-Cluster.
"""
const INPUTCLUSTER = parsed_args["inputCluster"] 

##################################################

include("Constants.jl") # Physical constants
include("Bath.jl") # Definition of the background bath
include("Mean.jl") # Mean quantities of the model
include("Cluster.jl") # Definition of the cluster stars
include("IntTable.jl") # Definition of the structure used for the computation of integrals K, dK
include("Klnnp.jl") # Computation of the Klnnp coefficients
include("ResonanceLines.jl") # Construction of the resonant lines
include("SRR.jl") # Computation of the SRR diffusion coefficients
include("NR.jl") # Computation of the NR diffusion coefficients
include("DJ.jl") # Computation of diffusion coefficient DRR_j (and DRR_jj in parallel)