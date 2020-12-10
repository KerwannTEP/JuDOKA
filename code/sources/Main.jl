include("Packages.jl") #import the packages
##################################################
include("Args.jl") # Parsing the command-line
##################################################
# General parameters
##################################################
"""
    PARALLEL
    
Determining if the code is run in parallel.
Set as `no` by default.
"""
const PARALLEL = parsed_args["parallel"]
if ((PARALLEL != "yes") && (PARALLEL != "no"))
    error("ERROR: UNKNOWN PARALLEL") # Unknown parallel procedure
end

"""
    lstart
    
Lower harmonic number used to compute `alpha(a,j)` the powerlaw index of the remainder of DRRjj.
Set at 1 by default.
"""
const lstart = parsed_args["lstart"]

"""
    lmax
    
Maximum harmonic number used in the multipole expansion of DRRjj.
Set as 10 by default.
"""
const lmax = parsed_args["lmax"]

"""
    nbK
    
Number of sampling points for geometric factor Klnnp.
Set as 100 by default.
"""
const nbK = parsed_args["nbK"]

"""
    nbResPoints
    
Number of points on the resonance lines.
Set as 100 by default.
"""
const nbResPoints = parsed_args["nbResPoints"]

"""
    mBH
    
Mass of the central BH, in solar masses (MSun).
Set as 4300000.0 MSun, the mass of the SMBH Sgr A*.
"""
const mBH = parsed_args["mBH"]

"""
    rh
    
Influence radius of the BH, in milliparsecs (mpc).
We do not compute any resonance lines beyond that point.
Set as 2000.0 mpc.
"""
const rh = parsed_args["rh"]

"""
    INPUTBATH
    
Name of the input file for the bath.
Set as `sources/inputBath.txt` by default.
"""
const INPUTBATH = parsed_args["inputBath"] 

"""
    INPUTCLUSTER

Name of the input file for the S-Cluster.
Set as `sources/inputCluster.txt` by default.
"""
const INPUTCLUSTER = parsed_args["inputCluster"] 

##################################################
include("Constants.jl") # Physical constants
include("Bath.jl") # Definition of the background bath
include("Cluster.jl") # Definition of the cluster stars

include("IntegrationTable.jl") # Definition of the structure used for the computation of integrals K, dK

include("Mean.jl") # Mean quantities of the model
include("Klnnp.jl") # Computation of the Klnnp coefficients
include("ResonanceLines.jl") # Construction of the resonant lines
include("SRR.jl") # Computation of the SRR diffusion coefficients
include("NR.jl") # Computation of the NR diffusion coefficients
include("DJ.jl") # Computation of diffusion coefficient DRR_j (and DRR_jj in parallel)