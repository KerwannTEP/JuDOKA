"""
    get_dataBath()
    
Returns the input data from the bath as an array of shape {nbBath,4}.
This array is stored in the order {gamma,m,a0,Mbar=M(<a0)}.
It is later used to get nbBath and to fill in tabgammaBar,tabmBath,taba0Bath,tabMbarBath,tabNbarBath.
"""
function get_dataBath()
    dataBath, headerBath = readdlm(INPUTBATH,header=true) # Reading the data file with the header
    return dataBath # Output the entire data
end
##################################################

"""
    dataBath
    
Input data from the Bath.
Use `get_dataBath()` to fill it in.
"""
const dataBath = get_dataBath()

"""
    nbBath
    
Number of bath populations.
Initialize `dataBath` with `get_dataBath()` first before looking at this value.
"""
const nbBath = size(dataBath)[1]

"""
    gBath(gamma)
    
Dimensionless function g(gamma), where gamma is the cusp of the family.
Used to make the translation between N(a) in SMA-space and N(<a) in radius space.
"""
function gBath(gamma::Float64)
    return 2.0^(-gamma)*sqrt(pi)*(SpecialFunctions.gamma(1.0+gamma))/(SpecialFunctions.gamma(gamma-0.5)) # Output of g(gamma)
end
##################################################
# Function that returns the normalisation prefactor
# Nbar in SMA-space, for a given set {gamma,m,a0,Mbar}
# for a given set {gamma,m,a0,Mbar}
##################################################
"""
    get_NbarBath(gamma,m,a0,Mbar)
    
Returns the normalization prefactor NBar in SMA-space, for a given set {gamma,m,a0,Mbar}.
The prefactor Nbar is such that `N(a)=Nbar*(a/a0)^(2-gamma)`.

# Arguments
- `gamma::Float64`: Cusp of the family.
- `m    ::Float64`: Individual masses of the family.
- `a0   ::Float64`: Scale radius of the family.
- `Mbar ::Float64`: Enclosed mass within a radius a0 of the family.
"""
function get_NbarBath(gamma::Float64,m::Float64,a0::Float64,Mbar::Float64)
    return (3.0-gamma)*gBath(gamma)*Mbar/(a0*m) # Value of Nbar such that N(a)=Nbar*(a/a0)^(2-gamma) in SMA-space
end
##################################################
"""
    tabgammaBath
    
Array containing the cusp index for all the families.
"""
const tabgammaBath = MVector{nbBath,Float64}([dataBath[iBath,1] for iBath=1:nbBath])

"""
    tabmBath
    
Array containing the individual masses for all the families.
"""
const tabmBath = MVector{nbBath,Float64}([dataBath[iBath,2] for iBath=1:nbBath])

"""
    taba0Bath
    
Array containing the scale radius `a0` for all the families.
"""
const taba0Bath = MVector{nbBath,Float64}([dataBath[iBath,3] for iBath=1:nbBath])

"""
    tabMbarBath
    
Array containing the enclosed masses, `Mbar=M(<a0)` for all the families.
"""
const tabMbarBath = MVector{nbBath,Float64}([dataBath[iBath,4] for iBath=1:nbBath])

"""
    tabNbarBath
    
Array containing the prefactors `Nbar` for all the families.
"""
const tabNbarBath = MVector{nbBath,Float64}([get_NbarBath(tabgammaBath[iBath],tabmBath[iBath],taba0Bath[iBath],tabMbarBath[iBath]) for iBath=1:nbBath]) 

##################################################

"""
    MenclosedBath(iBath,a)
    
Returns the enclosed mass `M(<a)` for the bath population indexed by iBath.
Equal to `M(<a)=Mbar*(a/a0)^(3-gamma)`.
"""
function MenclosedBath(iBath::Int64,a::Float64)
    gamma = tabgammaBath[iBath] # Cusp index of the current bath population
    a0 = taba0Bath[iBath] # Scale radius of the current bath population
    Mbar = tabMbarBath[iBath] # Enclosed mass M(<a0) for the current bath population
    return Mbar*(a/a0)^(3.0-gamma) # Output
end

"""
    NBath(iBath,a)
    
Returns the distribution function in SMA-space of `N(a)` for the bath component indexed by iBath.
Equal to `N(a)=Nbar*(a/a0)^(2-gamma)`.
"""
function NBath(iBath::Int64,a::Float64)
    gamma = tabgammaBath[iBath] # Cusp index of the current bath population
    Nbar  = tabNbarBath[iBath]
    a0    = taba0Bath[iBath] # Scale radius of the current bath population
    # Normalisation in SMA-space of the current population
    #####
    return Nbar*(a/a0)^(2.0-gamma) # Output
end

"""
    NBathTot(a)
    
Returns the total distribution function in SMA-space of `N(a)`.
"""
function NBathTot(a::Float64)
    res = 0.0
    for iBath=1:nbBath # Loop over the bath components
        m = tabmBath[iBath]^(2) # Individual mass^2 of the current bath component
        res += m*NBath(iBath,a) # Contribution from the bath component
    end
    return res
end

"""
    fjBath(a,j)
    
Returns the conditional PDF fj(j|a).
Normalized so that integrating `fj(j|a)` over `[jlc(a),1]` yields 1.
Proportional to `2j` in the range `[jlc(a),1]`.

# Remark:
- Independent of the family of the bath.
- Should add a test to ensure that one has `jlc<j<1.0`.
"""
function fjBath(a::Float64,j::Float64)
    return (2.0*j)/(1.0-(jlc(a))^(2)) # Output
end

"""
    dfjBathdj(a,j)
    
Returns the conditional PDF fj(j|a).
Normalized so that integrating `fj(j|a)` over `[jlc(a),1]` yields 1.

# Remark:
- Independent of the family of the bath.
- Should add a test to ensure that one has `jlc<j<1.0`.
"""
function dfjBathdj(a::Float64,j::Float64)
    return (2.0)/(1.0-(jlc(a))^(2))
end

"""
    FtotBath(a,j)

Wrapped function of the distributions functions which appears as an integrand in the DRR coefficients DRRJJ and DRRJ.
"""
function FtotBath(a::Float64,j::Float64)
    res = 0.0 # Initialising the result
    #####
    for iBath=1:nbBath # Loop over the bath components
        m = tabmBath[iBath] # Individual mass of the current bath component
        res += m^(2)*NBath(iBath,a) # Contribution from the bath component
    end
    res *= fjBath(a,j)
    #####
    return res # Output
end
