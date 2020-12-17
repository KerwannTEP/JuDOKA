##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--parallel"
    help = "Parallel computation: yes/no"
    arg_type = String
    default = "no"
    "--lmax"
    help = "Maximum harmonic number"
    arg_type = Int64
    default = 10
    "--nbK"
    help = "Number of sampling points for Klnnp"
    arg_type = Int64
    default = 100
    "--nbResPoints"
    help = "Number of points on the resonance lines"
    arg_type = Int64
    default = 100
    "--mBH"
    help = "Mass of the central BH (in Msun)"
    arg_type = Float64
    default = 4300000.0
    "--rh"
    help = "Influence radius of the BH (in mpc)"
    arg_type = Float64
    default = 2000.0
    "--inputBath"
    help = "Name of the input file for the bath"
    arg_type = String
    default = "../background/inputBath.txt"
    "--inputCluster"
    help = "Name of the input file for the S-cluster"
    arg_type = String
    default = "../background/inputCluster.txt"
end
parsed_args = parse_args(tabargs)
