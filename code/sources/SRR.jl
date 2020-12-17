"""
    DRR_JJnnp(n,np,a,j,IntTable_K=IntTable_K_serial,tabSMARes=tabSMARes_serial)

Computation of DRR_JJnn, contribution of the (n,np) resonance number to DRR.

Using log-sampling in the integration for better interpolation.

# Arguments
- `n   ::Int64`  : first resonance number of the orbit.
- `np  ::Int64`  : second resonance number of the orbit.
- `a   ::Float64`: semi-major axis of the orbit.
- `j   ::Float64`: reduced angular momentum of the orbit.
- `IntTable_K=IntTable_K_serial:: IntTable`: Structure which contains the tables needed for the computation of this integral.
- `tabSMARes=tabSMARes_serial`: Static container for the resonant lines, contains the edge of the domain [amin,amax] of the resonant line.
"""
function DRR_JJnnp(n::Int64,np::Int64,a::Float64,j::Float64,IntTable_K::IntTable=IntTable_K_serial,tabSMARes::Array{Float64,2}=tabSMARes_serial)
    #####
    omega = n*nup(a,j) # Value of the resonant frequency
    #####
    ResPoints = tabSMARes!(n,np,a,j,tabSMARes) # Trying to compute the resonant line
    #####
    if (ResPoints == 0)
        return 0.0 # If the resonant line is empty, we return 0.0 for the diffusion coefficient
    end
    #####
    # First resonance line
    apmin1, apmax1 = tabSMARes[1,1], tabSMARes[1,2] # Boundary of the SMA-region where the resonant line lies
    dalphap1 = log(apmax1/apmin1)/(nbResPoints) # Step distance in log-space
    
    if (ResPoints == 2) # If second resonance line
        apmin2, apmax2 = tabSMARes[2,1], tabSMARes[2,2] # Boundary of the SMA-region where the resonant line lies
        dalphap2 = log(apmax2/apmin2)/(nbResPoints) # Step distance in log-space
    end
    #####
    res = 0.0 # Initialisation of the result
    #####
    for alphap1 = range(log(apmin1),length=nbResPoints,log(apmax1)-dalphap1) # Uniform logarithmic spacing of the SMA !! ATTENTION, we do not go through the last element, as we will sum in the middle of the interval
        ap = exp(alphap1 + 0.5*dalphap1) # Central position of the SMA considered
        jp = getjRes(np,ap,omega) # Finding the location of the resonant jp
        res += dalphap1*ap*FtotBath(ap,jp)*AnnpSQ(n,np,a,j,ap,jp,IntTable_K)/(abs(dnupdj(ap,jp))) # Local contribution !! ATTENTION, not to forget the ap
    end
    if (ResPoints == 2) # If two resonance lines
        for alphap2 = range(log(apmin2),length=nbResPoints,log(apmax2)-dalphap2) # Uniform logarithmic spacing of the SMA !! ATTENTION, we do not go through the last element, as we will sum in the middle of the interval
            ap = exp(alphap2 + 0.5*dalphap2) # Central position of the SMA considered
            jp = getjRes(np,ap,omega) # Finding the location of the resonant jp
            res += dalphap2*ap*FtotBath(ap,jp)*AnnpSQ(n,np,a,j,ap,jp,IntTable_K)/(abs(dnupdj(ap,jp))) # Local contribution !! ATTENTION, not to forget the ap
        end
    end

    #####
    res *= (n^(2))/(abs(np)) # Adding the prefactor in n^2/|np|
    return res # Output
end

"""
    DRR_JJ(a,j,IntTable_K=IntTable_K_serial,tabSMARes=tabSMARes_serial)

Computation of DRR_JJ the scalar resonant diffusion coefficient in JJ-coordinates.

# Arguments
- `a   ::Float64`: semi-major axis of the orbit.
- `j   ::Float64`: reduced angular momentum of the orbit.
- `IntTable_K=IntTable_K_serial:: IntTable`: Structure which contains the tables needed for the computation of this integral.
- `tabSMARes=tabSMARes_serial`: Static container for the resonant lines, contains the edge of the domain [amin,amax] of the resonant line.

"""
function DRR_JJ(a::Float64,j::Float64,IntTable_K::IntTable=IntTable_K_serial,tabSMARes::Array{Float64,2}=tabSMARes_serial)
    res = 0.0 # Initialisation of the result
    #####
    for iPair=1:nbResPairs # Loop over the resonant pairs
        n, np = tabResPairs[1,iPair], tabResPairs[2,iPair] # Value of the considered resonance (n,np)
        res += DRR_JJnnp(n,np,a,j,IntTable_K,tabSMARes) # Adding the contribution from the resonance pair (n,np)
    end
    #####
    res *= 4.0*pi*G^(2) # Adding the prefactor 4piG^2
    return res # Output
end

"""
    DRR_jj(a,j,IntTable_K=IntTable_K_serial,tabSMARes=tabSMARes_serial)

Computation of DRR_jj the (normalized) scalar resonant diffusion coefficient in jj-coordinates.

# Arguments
- `a   ::Float64`: semi-major axis of the orbit.
- `j   ::Float64`: reduced angular momentum of the orbit.
- `IntTable_K=IntTable_K_serial:: IntTable`: Structure which contains the tables needed for the computation of this integral.
- `tabSMARes=tabSMARes_serial`: Static container for the resonant lines, contains the edge of the domain [amin,amax] of the resonant line.

"""
function DRR_jj(a::Float64,j::Float64,IntTable_K::IntTable=IntTable_K_serial,tabSMARes::Array{Float64,2}=tabSMARes_serial)
    return DRR_JJ(a,j,IntTable_K,tabSMARes)/((Jc(a))^(2)) # Output
end