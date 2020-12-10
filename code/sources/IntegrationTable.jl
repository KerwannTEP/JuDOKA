##################################################
# Structure which gathers the quantities used in computing K and dK
##################################################
"""
Structure which gathers all the necessary arrays needed to compute the geometrical factors (Knn and its derivatives).

# Summary

struct IntTable <: Any

# Fields
- `tabR      :: MArray{Tuple{nbK},Float64,1,nbK}`  : Table for the radius of the first particle.
- `tabg      :: MArray{Tuple{nbK},Float64,1,nbK}`  : Table for the prefactor of the first particle.
- `tabRp     :: MArray{Tuple{nbK},Float64,1,nbK}`  : Table for the radius of the second particle.
- `tabgp     :: MArray{Tuple{nbK},Float64,1,nbK}`  : Table for the prefactor of the second particle.
- `tabw      :: MArray{Tuple{nbK+2},Int64,1,nbK+2}`: Table for the sorting of the arrays.
- `tabh1     :: MArray{Tuple{nbK},Float64,1,nbK}`  : Table for the h1 of the first particle.
- `tabh2     :: MArray{Tuple{nbK},Float64,1,nbK}`  : Table for the h2 of the second particle.
- `tabdeltaP :: MArray{Tuple{nbK},Float64,1,nbK}`  : Table of the partial sum deltaP.
- `tabdeltaQ :: MArray{Tuple{nbK},Float64,1,nbK}`  : Table of the partial sum deltaQ.

"""
struct IntTable
    tabR::MArray{Tuple{nbK},Float64,1,nbK}	# Table for the radius of the first particle
    tabg::MArray{Tuple{nbK},Float64,1,nbK}	# Table for the prefactor of the first particle
    tabRp::MArray{Tuple{nbK},Float64,1,nbK}	# Table for the radius of the second particle
    tabgp::MArray{Tuple{nbK},Float64,1,nbK}	# Table for the prefactor of the second particle
    tabw::MArray{Tuple{nbK+2},Int64,1,nbK+2}	# Table for the sorting of the arrays
    tabh1::MArray{Tuple{nbK},Float64,1,nbK}	# Table for the h1 of the first particle
    tabh2::MArray{Tuple{nbK},Float64,1,nbK}	# Table for the h2 of the second particle
    tabdeltaP::MArray{Tuple{nbK},Float64,1,nbK}	# Table of the partial sum deltaP
    tabdeltaQ::MArray{Tuple{nbK},Float64,1,nbK} # Table of the partial sum deltaQ
end
##################################################
# Creation of a structure
##################################################

"""
    tab1

Vector of zeroes of length `nbK` used in the initialization of `IntTable` structures.
"""
const tab1 = @MVector zeros(Float64,nbK)

"""
    tab2

Vector of zeroes of length `nbK+2` used in the initialization of `IntTable` structures.
"""
const tab2 = @MVector zeros(Int64,nbK+2)

"""
    Table_create!(t1=tab1, t2=tab2)

Creates a IntTable structure with all fields set as zeros-arrays.
"""
function Table_create!(t1=tab1, t2=tab2)
    return IntTable(copy(t1),copy(t1),copy(t1),copy(t1),copy(t2),copy(t1),copy(t1),copy(t1),copy(t1))
end
##################################################
# Initialization of a structure
##################################################

"""
    Table_init!(IntTable_arg, <keyword arguments>)

Initializes the IntTable structure `IntTable_arg` at the given arguments.

# Arguments
- `n  ::Int64`  : resonance number of the first orbit.
- `np ::Int64`  : resonance number of the second orbit.
- `a  ::Float64`: semi-major axis of the first orbit.
- `j  ::Float64`: reduced angular momentum of the first orbit.
- `ap ::Float64`: semi-major axis of the decond orbit.
- `jp ::Float64`: reduced angular momentum of the second orbit.

"""
function Table_init!(IntTable_arg,n::Int64,np::Int64,a::Float64,j::Float64,ap::Float64,jp::Float64)
    for k=1:nbK # Loop over the sample points
        cosf = tabcosf_Klnnp[k] # Extracting the value of cos(f)
        cosnf, cosnpf = tabcosnf_Klnnp[abs(n),k], tabcosnf_Klnnp[abs(np),k] # Extracting the value of cos(n*f) and cos(np*f) !! ATTENTION, we use (abs(n),abs(np)) for the array access to be ok.
        R, Rp = a*(j^(2))/(1.0+sqrt(1.0-j^(2))*cosf), ap*(jp^(2))/(1.0+sqrt(1.0-jp^(2))*cosf) # Radii of the first and second particles
        dMdf, dMpdfp = (R/a)^(2)/(j), (Rp/ap)^(2)/(jp) # Computing the Jacobian dM/df and dMp/dfp
        dM = ddMdfdj(j,cosf) # Computing j-derivative of the Jacobian
        ps = psi(j,cosf) # Computing psi(j,f)
        dps = dpsidj(j,cosf) # Computing j-derivative of psi(j,f)
        #####
	IntTable_arg.tabR[k], IntTable_arg.tabg[k]= R, dMdf*cosnf
	IntTable_arg.tabRp[k], IntTable_arg.tabgp[k] = Rp, dMpdfp*cosnpf
	IntTable_arg.tabh1[k], IntTable_arg.tabh2[k] = dM*cosnf, ps*dps/j *cosnf #h1, h2
    end
    IntTable_arg.tabw[1] = 0 # Initial term. !! ATTENTION, indices are shifted by one
    icount = 0 # Temporary counter for filling on tabw
    for k=1:nbK
	while (icount < nbK && IntTable_arg.tabR[icount+1] <= IntTable_arg.tabRp[k]) # !! ATTENTION, to the order of the tests to avoid segmentation faults
            icount += 1 # Updating the counter
        end
        IntTable_arg.tabw[k+1] = icount # ATTENTION, indices are shifted by one
    end
    IntTable_arg.tabw[nbK+2] = nbK # Last term. !! ATTENTION, indices are shifted by one
end
