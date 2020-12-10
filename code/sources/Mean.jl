"""
    tabylnSQ
    
Container for the squared spherical harmonics prefactor |Y_l^n|^2.
ATTENTION, harmonic index starts at l = 0.
We pre-compute all the values of |Y_l^n (pi/2,pi/2)|^2 rather having to load any special functions librairies in order to save time.

Use `tabylnSQ!()` to fill it in.
"""
const tabylnSQ = zeros(Float64,lmax+1,lmax+1)

"""
    tabylnSQ!()
    
Fills the array (of size lmax*lmax) tabylnSQ with the values of the |Y_l^n|^2, the squared spherical harmonics prefactors.
tabylnSQ is indexed by two integers [l,n], with l,n going from 1 to lmax+1 since arrays begin at 1.
tabylnSQ[l+1,n+1] contains the value `|Y_l^n|^2 = |Yln(pi/2,pi/2)|^2`. 
Those are computed by recursion.
"""
function tabylnSQ!()
    # Making sure that the array is initially filled with zeros
    for n=1:(lmax+1), l=1:(lmax+1) # ATTENTION, the size of the array is (lmax+1) in both directions
        tabylnSQ[l,n] = 0.0
    end
    #####
    # Computing all the squared elements by recurrence !! ATTENTION, this corresponds to the squared elements
    tabylnSQ[1,1] = 1.0/(4.0*pi) # |Yln(pi/2,pi/2)|^(2) for (l,n)=(0,0) !! HARMONIC INDEX STARTS AT \ell=0
    # First loop to fill in the diagonal elements
    for l=0:(lmax-1) #!! ATTENTION TO THE BOUNDS OF THE LOOP
        tabylnSQ[(l+1)+1,(l+1)+1] = ((2.0*l+3.0)/(2.0*l+2.0))*tabylnSQ[l+1,l+1]
    end
    # Second loop to fill in the remaining terms
    for n=0:lmax # We perform the loop at fixed m
        for l=n:2:(lmax-2) # At fixed n, we scan the \ell by step of 2 !! ATTENTION TO THE BOUNDS OF THE LOOP
            tabylnSQ[(l+2)+1,n+1] = ((2.0*l+5.0)/(2.0*l+1.0))*(((l-n+1.0)*(l+n+1.0))/((l-n+2.0)*(l+n+2.0)))*tabylnSQ[l+1,n+1]
        end
    end
end
##################################################
tabylnSQ!() # Initialising tabylnSQ with the pre-computed coefficients
##################################################
"""
    ylnSQ(l,n)

Wrapped function that returns the squared spherical harmonics prefactor |Y_l^n (pi/2,0)|^2.
Caution: this is the SQUARED quantity.
"""
function ylnSQ(l::Int64,n::Int64)
    return tabylnSQ[l+1,abs(n)+1] # Using the symmetry |Yl,n(pi/2,pi/2)|^2 = |Yl,-n(pi/2,pi/2)|^2
end

##################################################
# Definitions of useful quantities
##################################################

"""
    rg
    
Schwartzschild (gravitational) radius. 
Defined as `(G*mBH)/(cvel^(2))`, where `G` is Newton's constant, `mBH` the mass of the SMBH and `cvel` the light celerity.
"""
const rg = (G*mBH)/(cvel^(2))

"""
    Jc(a)
    
Angular momentum of a circular orbit with radius a.
"""
function Jc(a::Float64)
    return sqrt(G*mBH*a) 
end

"""
    jlc(a)
    
Reduced angular momentum of the loss-cone for a semi-major axis a.
For j<jlc(a), there are no allowed orbits
"""
function jlc(a::Float64)
    return 4.0*sqrt(rg/a) 
end

"""
    nuKep(a)
    
Fast Keplerian frequency of an orbit with semi-major axis a.
"""
function nuKep(a::Float64)
    return sqrt((G*mBH)/(a^(3))) 
end

"""
    nuGR(a,j)
    
Schwarzschild precession frequency for an orbit of semi-major axis a and reduced angular momentum J/Jc = j.
"""
function nuGR(a::Float64,j::Float64)
    return 3.0*(rg/a)*nuKep(a)/(j^(2)) 
end

"""
    dnuGRdj(a,j)
    
j-partial gradient of the Schwarzschild precession frequency for an orbit of semi-major axis a and reduced angular momentum J/Jc = j.
"""
function dnuGRdj(a::Float64,j::Float64)
    return -6.0*(rg/a)*nuKep(a)/(j^(3))
end

"""
    d2nuGRdj2(a,j)
    
jj-second-order gradient of the Schwarzschild precession frequency for an orbit of semi-major axis a and reduced angular momentum J/Jc = j.
"""
function d2nuGRdj2(a::Float64,j::Float64)
    return 18.0*(rg/a)*nuKep(a)/(j^(4))
end

"""
    LegendreP(n,x)
   
Wrapped function to evalue the Legendre functions with non-integer n at argument x.
This uses the Hypergeometric function expression `Pn(x) = ₂F₁(-n,n+1.0,1.0,(1.0-x)*0.5) `.
"""
function LegendreP(n::Float64,x::Float64)
    return _₂F₁(-n,n+1.0,1.0,(1.0-x)*0.5) 
end

##################################################
# Computation of actual physically relevant functions
##################################################

"""
    nuMbar(M_below_a,a)
    
SMA-dependent part of themass precession frequency.
Defined as `nuKep*M(<a)/mBH`, where `M(<a)` is the mass enclosed within semi-major axis a and `mBH` the mass of the SMBH.
"""
function nuMbar(M_below_a::Float64,a::Float64)
    return nuKep(a)*(M_below_a)/(mBH) # Output
end

##################################################

"""
    nbjInt
    
Number of j used for the interpolation of the functions j->hM(j) and its derivatives.
Set as 1000 by default.
"""
const nbjInt = 1000

"""
    jmaxhM
    
For j>jmaxhM, we use a linear approximation for j->hM(j) and its derivatives.
Set as 0.99 by default.
"""
const jmaxhM = 0.99

##################################################

"""
    _hM(gamma,j,p1,p2)
    
j-dependent part of the mass precession frequency.
We should take p1 = LegendreP[1-gamma,1/j] and p2 = LegendreP[2-gamma,1/j].
Depends on the cusp index gamma and the reduced angular momentum j.

# Remarks
- We use a quadratic approximation for j>jmaxhM 
- This function is not continuous since we use a linear approximation for j>jmaxhM. By default, jmaxhM=0.99.
"""
function _hM(gamma::Float64,j::Float64,p1::Float64,p2::Float64)
    if (j > jmaxhM) # We use the linear approximation
        return (1.0/2.0)*(-3.0+gamma)-
               (1.0/8.0)*(-12.0+gamma+4.0*gamma^(2)-gamma^(3))*(1.0-j)+
               (1.0/96.0)*(-3.0+gamma)*(-1.0+gamma)*gamma*(-16.0+gamma*(3.0+gamma))*(1.0-j)^(2) # Output
    else # We can use the generic expression
        return (j^(4.0-gamma))/(1.0-j^(2))*(p1 - p2/(j)) # Output
    end
end

"""
    _dhMdj(gamma,j,p1,p2)
    
j-gradient of the j-dependent part of the mass precession frequency.
We should take p1 = LegendreP[1-gamma,1/j] and p2 = LegendreP[2-gamma,1/j].
Depends on the cusp index gamma and the reduced angular momentum j.

# Remarks
- We use a quadratic approximation for j>jmaxhM 
- This function is not continuous since we use a linear approximation for j>jmaxhM. By default, jmaxhM=0.99.
"""
function _dhMdj(gamma::Float64,j::Float64,p1::Float64,p2::Float64)
    if (j > jmaxhM) # We use the linear approximation
        return (1.0/8.0)*(-12.0+gamma+4.0*gamma^(2)-gamma^(3))-
               (1.0/48.0)*(-3.0+gamma)*(-1.0+gamma)*gamma*(-16.0+gamma*(3.0+gamma))*(1.0-j)-
               (1.0/768)*((-3.0+gamma)*(-1.0+gamma)*gamma*(72.0+gamma*(-94+gamma*(-13.0+gamma*(10.0+gamma)))))*(1.0-j)^(2) # Output
    else # We can use the generic expression
        return (j^(2.0-gamma))/((1.0-j^(2))^(2))*(j*(2.0+(2.0-gamma)*(1.0-j^(2)))*p1 - (1.0+j^(2))*p2) # Output
    end
end

"""
    _d2hMdj2(gamma,j,p1,p2)
    
Second-order j-gradient of the j-dependent part of the mass precession frequency.
We should take p1 = LegendreP[1-gamma,1/j] and p2 = LegendreP[2-gamma,1/j].
Depends on the cusp index gamma and the reduced angular momentum j.

# Remarks
- We use a quadratic approximation for j>jmaxhM 
- This function is not continuous since we use a linear approximation for j>jmaxhM. By default, jmaxhM=0.99.
"""
function _d2hMdj2(gamma::Float64,j::Float64,p1::Float64,p2::Float64)
    if (j > jmaxhM) # We use the linear approximation
        return (1.0/48.0)*(-3.0+gamma)*(-1.0+gamma)*gamma*(-16.0+gamma*(3.0+gamma))+
               (1.0/384.0)*(-3.0+gamma)*(-1.0+gamma)*gamma*(72.0+gamma*(-94.0+gamma*(-13.0+gamma*(10.0+gamma))))*(1.0-j)+
               (1.0/7680.0)*(-3.0+gamma)*(-1.0+gamma)*gamma*(1.0+gamma)*(2.0+gamma)*(360.0+gamma*(-358.0+gamma*(19.0+gamma*(18.0+gamma))))*(1.0-j)^(2) # Output
    else # We can use the generic expression
        return (j^(2.0-gamma))/((1.0-j^(2))^(3))*((2.0*(3.0-gamma)^(2)-3.0*j^(2)*(4.0-gamma)*(1.0-gamma)+j^(4)*(2.0-gamma)*(1.0-gamma))*p1+
                                                  j*(-12.0+j^(2)*(4.0-gamma)*(1.0-gamma)+(5.0-gamma)*gamma)*p2) # Output
    end
end

"""
    gethMInt(iBath)
    
Function that returns an interpolation function for the function j->hM(j) for the bath population indexed by iBath in the source file.
hM(j) is the j-dependent part of the mass precession frequency.
"""
function gethMInt(iBath::Int64)
    jminInt = jlc(rh) # Minimum value of j for the interpolation, fixed by the influence radius
    jmaxInt = 1.0 # Maximum value of j for the interpolation
    #####
    rangejInt = range(jminInt,length=nbjInt,jmaxInt) # Range in j over which the interpolation is performed
    tabjInt = collect(rangejInt) # Table of the values of j for which the functions must be pre-computed
    tabhMInt = zeros(Float64,nbjInt) # Container for the values of hM(j)
    ####
    gamma = tabgammaBath[iBath] # Cusp index of the current bath population
    #####
    for indj=1:nbjInt # Loop over the j of the grid
        jloc = tabjInt[indj] # Current value of j
        p1 = LegendreP(1.0-gamma,1.0/(jloc)) # Current value of LegendreP[1-gamma,1/j]
        p2 = LegendreP(2.0-gamma,1.0/(jloc)) # Current value of LegendreP[2-gamma,1/j]
        tabhMInt[indj] = _hM(gamma,jloc,p1,p2) # Pre-computing the value of hM(j)
    end
    #####
    inthM = Interpolations.scale(interpolate(tabhMInt, BSpline(Cubic(Line(OnGrid())))),rangejInt) # Constructing the interpolation function for j->hM(j)
    #####
    return inthM # Returning the interpolation function
end

"""
    getdhMdjInt(iBath)
    
Function that returns an interpolation function for the function j->dhMdj(j) for the bath population indexed by iBath in the source file.
dhMdj(j) is the j-gradient of the j-dependent part of the mass precession frequency.
"""
function getdhMdjInt(iBath::Int64)
    jminInt = jlc(rh) # Minimum value of j for the interpolation, fixed by the influence radius
    jmaxInt = 1.0 # Maximum value of j for the interpolation
    #####
    rangejInt = range(jminInt,length=nbjInt,jmaxInt) # Range in j over which the interpolation is performed
    tabjInt = collect(rangejInt) # Table of the values of j for which the functions must be pre-computed
    tabdhMdjInt = zeros(Float64,nbjInt) # Container for the values of dhMdj(j)
    ####
    gamma = tabgammaBath[iBath] # Cusp index of the current bath population
    #####
    for indj=1:nbjInt # Loop over the j of the grid
        jloc = tabjInt[indj] # Current value of j
        p1 = LegendreP(1.0-gamma,1.0/(jloc)) # Current value of LegendreP[1-gamma,1/j]
        p2 = LegendreP(2.0-gamma,1.0/(jloc)) # Current value of LegendreP[2-gamma,1/j]
        tabdhMdjInt[indj] = _dhMdj(gamma,jloc,p1,p2) # Pre-computing the value of hM(j)
    end
    #####
    intdhMdj = Interpolations.scale(interpolate(tabdhMdjInt, BSpline(Cubic(Line(OnGrid())))),rangejInt) # Constructing the interpolation function for j->hM(j)
    #####
    return intdhMdj # Returning the interpolation function
end

"""
    getd2hMdj2Int(iBath)
    
Function that returns an interpolation function for the function j->d2hMdj2(j) for the bath population indexed by iBath in the source file.
d2hMdj2(j) is the jj-second-order gradient of the j-dependent part of the mass precession frequency.
"""
function getd2hMdj2Int(iBath::Int64)
    jminInt = jlc(rh) # Minimum value of j for the interpolation, fixed by the influence radius
    jmaxInt = 1.0 # Maximum value of j for the interpolation
    #####
    rangejInt = range(jminInt,length=nbjInt,jmaxInt) # Range in j over which the interpolation is performed
    tabjInt = collect(rangejInt) # Table of the values of j for which the functions must be pre-computed
    tabd2hMdj2Int = zeros(Float64,nbjInt) # Container for the values of dhMdj(j)
    ####
    gamma = tabgammaBath[iBath] # Cusp index of the current bath population
    #####
    for indj=1:nbjInt # Loop over the j of the grid
        jloc = tabjInt[indj] # Current value of j
        p1 = LegendreP(1.0-gamma,1.0/(jloc)) # Current value of LegendreP[1-gamma,1/j]
        p2 = LegendreP(2.0-gamma,1.0/(jloc)) # Current value of LegendreP[2-gamma,1/j]
        tabd2hMdj2Int[indj] = _d2hMdj2(gamma,jloc,p1,p2) # Pre-computing the value of hM(j)
    end
    #####
    intd2hMdj2 = Interpolations.scale(interpolate(tabd2hMdj2Int, BSpline(Cubic(Line(OnGrid())))),rangejInt) # Constructing the interpolation function for j->hM(j)
    #####
    return intd2hMdj2 # Returning the interpolation function
end

##################################################
# Preparing containers for the interpolation function of j->hM(j) and j->dhMdj(j)
# See https://discourse.julialang.org/t/array-of-functions-is-there-a-way-to-avoid-allocations-performance-penalty/24471

"""
    tabhMInt
    
Static container for the interpolation function of j->hM(j).
Set equal to `SVector{nbBath}([gethMInt(iBath) for iBath=1:nbBath])`.
See `https://discourse.julialang.org/t/array-of-functions-is-there-a-way-to-avoid-allocations-performance-penalty/24471`.
"""
const tabhMInt = SVector{nbBath}([gethMInt(iBath) for iBath=1:nbBath])

"""
    tabdhMIntdj
    
Static container for the interpolation function of j->dhMdj(j).
Set equal to `SVector{nbBath}([getdhMdjInt(iBath) for iBath=1:nbBath])`.
See `https://discourse.julialang.org/t/array-of-functions-is-there-a-way-to-avoid-allocations-performance-penalty/24471`.
"""
const tabdhMdjInt = SVector{nbBath}([getdhMdjInt(iBath) for iBath=1:nbBath])

"""
    tabd2hMIntdj2
    
Static container for the interpolation function of j->d2hMdj2(j).
Set equal to `SVector{nbBath}([getd2hMdj2Int(iBath) for iBath=1:nbBath])`.
See `https://discourse.julialang.org/t/array-of-functions-is-there-a-way-to-avoid-allocations-performance-penalty/24471`.
"""
const tabd2hMdj2Int = SVector{nbBath}([getd2hMdj2Int(iBath) for iBath=1:nbBath])

##################################################

"""
    nuMBath(iBath,a,j)
    
Mass precession frequency for an orbit with semi-major axis a and reduced angular momentum j, for the family indexed by iBath in the source file.
Uses pre-computed interpolations functions stacked into arrays, as well as an @inline function for the evaluation to be fast.
"""
@inline function nuMBath(iBath::Int64,a::Float64,j::Float64)
    Menclosed = MenclosedBath(iBath,a) # Enclosed mass M(<a) of the current bath population
    #####
    return nuMbar(Menclosed,a)*tabhMInt[iBath](j) # Returning the mass precession frequency !! ATTENTION, we use the interpolation function for j->hM(j)
end

"""
    dnuMBathdj(iBath,a,j)
    
j-gradient of the mass precession frequency for an orbit with semi-major axis a and reduced angular momentum j, for the family indexed by iBath in the source file.
Uses pre-computed interpolations functions stacked into arrays, as well as an @inline function for the evaluation to be fast.
"""
@inline function dnuMBathdj(iBath::Int64,a::Float64,j::Float64)
    Menclosed = MenclosedBath(iBath,a) # Enclosed mass M(<a) of the current bath population
    #####
    return nuMbar(Menclosed,a)*tabdhMdjInt[iBath](j) # Returning the mass precession frequency !! ATTENTION, we use the interpolation function for j->dhMdj(j)
end

"""
    d2nuMBathdj2(iBath,a,j)
    
jj-second-ordre gradient of the mass precession frequency for an orbit with semi-major axis a and reduced angular momentum j, for the family indexed by iBath in the source file.
Uses pre-computed interpolations functions stacked into arrays, as well as an @inline function for the evaluation to be fast.
"""
@inline function d2nuMBathdj2(iBath::Int64,a::Float64,j::Float64)
    Menclosed = MenclosedBath(iBath,a) # Enclosed mass M(<a) of the current bath population
    #####
    return nuMbar(Menclosed,a)*tabd2hMdj2Int[iBath](j) # Returning the mass precession frequency !! ATTENTION, we use the interpolation function for j->dhMdj(j)
end

"""
    nuMtot(a,j)
    
Total mass precession frequency for an orbit with semi-major axis a and reduced angular momentum j.
Uses pre-computed interpolations functions stacked into arrays, as well as an @inline function for the evaluation to be fast.
"""
@inline function nuMtot(a::Float64,j::Float64)
    res = 0.0 # Initialising the result
    #####
    for iBath=1:nbBath # Loop over the bath particles
        res += nuMBath(iBath,a,j) # Contribution from the ith bath population
    end
    #####
    return res # Output
end

"""
    dnuMtotdj(a,j)
    
j-gradient of the total mass precession frequency for an orbit with semi-major axis a and reduced angular momentum j.
Uses pre-computed interpolations functions stacked into arrays, as well as an @inline function for the evaluation to be fast.
"""
@inline function dnuMtotdj(a::Float64,j::Float64)
    res = 0.0 # Initialising the result
    #####
    for iBath=1:nbBath # Loop over the bath particles
        res += dnuMBathdj(iBath,a,j) # Contribution from the ith bath population
    end
    #####
    return res # Output
end

"""
    d2nuMtotdj2(a,j)
    
jj-second-order gradient of the total mass precession frequency for an orbit with semi-major axis a and reduced angular momentum j.
Uses pre-computed interpolations functions stacked into arrays, as well as an @inline function for the evaluation to be fast.
"""
@inline function d2nuMtotdj2(a::Float64,j::Float64)
    res = 0.0 # Initialising the result
    #####
    for iBath=1:nbBath # Loop over the bath particles
        res += d2nuMBathdj2(iBath,a,j) # Contribution from the ith bath population
    end
    #####
    return res # Output
end

"""
    nup(a,j)
    
Total precession frequency for an orbit with semi-major axis a and reduced angular momentum j.
Uses pre-computed interpolations functions stacked into arrays, as well as an @inline function for the evaluation to be fast.
"""
@inline function nup(a::Float64,j::Float64)
    return nuGR(a,j) + nuMtot(a,j) # Returning the value of the total precession frequency
end

"""
    dnupdj(a,j)
    
j-gradient of the total precession frequency for an orbit with semi-major axis a and reduced angular momentum j.
Uses pre-computed interpolations functions stacked into arrays, as well as an @inline function for the evaluation to be fast.
"""
@inline function dnupdj(a::Float64,j::Float64)
    return dnuGRdj(a,j) + dnuMtotdj(a,j) # Returning the gradient of the total precession frequency
end

"""
    d2nupdj2(a,j)
    
jj-second-order gradient of the total precession frequency for an orbit with semi-major axis a and reduced angular momentum j.
Uses pre-computed interpolations functions stacked into arrays, as well as an @inline function for the evaluation to be fast.
"""
@inline function d2nupdj2(a::Float64,j::Float64)
    return d2nuGRdj2(a,j) + d2nuMtotdj2(a,j) # Returning the second gradient of the total precession frequency
end
