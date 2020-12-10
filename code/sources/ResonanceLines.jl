"""
    get_nbResPairs()
    
Returns the number of resonant pairs (n,np).
Here, the resonant pairs (n,np) must satisfy the 3 conditions:
- `l <= n <= lmax`.
- `-lmax <= np <= lmax`.
- `n` and `np` have the same parity.
"""
function get_nbResPairs()
    s = 0 # Initialising the counter
    for n=1:lmax # Loop over the index n !! ATTENTION, we have have n>0
        for np=-lmax:lmax # Loop over the index np
            if (np != 0 && mod(n-np,2) == 0) # For np to be allowed, we must have np != 0, and n,np of the same parity
                s += 1 # Updating the counter
            end
        end
    end
    return s # Output
end

##################################################

"""
    nbResPairs
    
Total number of resonant pairs (n,np).
Use `get_nbResPairs()` to determine this quantity.
"""
const nbResPairs = get_nbResPairs()

"""
    tabResPairs
    
Static container for the resonance pairs (n,np)
Use `tabResPairs!()` to fill it in.
"""
const tabResPairs = zeros(Int64,2,nbResPairs)

##################################################

"""
    tabResPairs!()
    
Fills in the table tabResPairs, of size 2*nbResPairs, of the resonance pairs (n,np) 
"""
function tabResPairs!()
    iPair = 1 # Counter for the resonant pair
    for n=1:lmax # Loop over the index n !! ATTENTION, we have have n>0
        for np=-lmax:lmax # Loop over the index np
            if (np != 0 && mod(n-np,2) == 0) # For np to be allowed, we must have np != 0, and n,np of the same parity
                tabResPairs[1,iPair], tabResPairs[2,iPair] = n, np # Filling in the array with (n,np)
                iPair += 1 # Updating the pair
            end
        end
    end
end

##################################################
tabResPairs!() # Preparing the table of the resonance pairs to consider
##################################################

"""
    bisection(fun, xl, xu, <keyword arguments>)

Algorithm which looks for zeroes of the functions `fun` given in argument.
It is a simple bisection algorithm, but it makes no allocations and is sufficiently fast.
It allows us not to have to use the Roots library that makes a lot of allocations.
The optional tolerances are set to the same as the ones found in Roots.jl, and are most likely overkill.

The zero is searched within the interval `[xl, xu]`, whose bounds are given as arguments.
"""
function bisection(fun, xl::Float64, xu::Float64, tolx::Float64=4.0*eps(Float64), tolf::Float64=4.0*eps(Float64), iterMAX::Int64=50)
    if (xl > xu)
        xl, xu = xu, xl # Ordering the input arguments
    end
    #####
    fl, fu = fun(xl), fun(xu) # Evaluating the function on the bracket
    #####
    @assert fl*fu < 0.0 "bisect: NOT A BRACKET"
    #####
    if (fl == 0.0) # We have already found a solution on the left bracket
        return xl # Returning the left bracket
    end
    #####
    if (fu == 0.0) # We have already found a solution on the right bracket
        return xu # Returning the right bracket
    end
    #####
    iter = 0 # Counter for the iterations
    #####
    while true # Bisection loop
        xm = (xl+xu)*0.5 # Middle value
        fm = fun(xm) # Value at the midpoint
        #####
        iter += 1 # Updating the iteration counter
        #####
        if (fm*fl < 0.0) # One root can be found between the left point and the middle point
            xu, fu = xm, fm # The upper point becomes the midpoint
        else
            xl, fl = xm, fm # The lower point becomes the midpoint
        end
        #####
        if ((abs(xu-xl) < tolx) || (abs(fm) <= tolf) || (iter > iterMAX)) # Testing whether or not we have found a root given the tolerances or if we have made enough iterations
            break # We have found a root
        end
    end
    return (xl+xu)*0.5 # Returning the central location of the bracket
end

##################################################

"""
    tabSMARes
    
Static container for the resonant lines, contains the edge of the domain [amin,amax] of the resonant line.
Use `tabSMARes!()` to fill it in.
"""
const tabSMARes = zeros(Float64,4)
"""
    amin
    
Minimum semi-major axis for which the resonance condition is solved, corresponds to jlc(amin)=0.97 !!
ATTENTION, we make sure that jlc(amin)<1.0 to always have meaningful brackets in SMAResQ().

Set as `17*rg` where `rg` is the Schwartzschild radius.
"""
const amin = 17.0*rg # Minimum SMA for which the resonance condition is solved, corresponds to jlc(amin)=0.97 !! ATTENTION, we make sure that jlc(amin)<1.0 to always have meaningful brackets in SMAResQ()

"""
    amax
    
Maximum semi-major axis for which the resonance condition is solved
Set as `rh` where `rh` is given as an argument of the initial command line as `--rh`.
"""
const amax = rh 

"""
    nbSMAResSEARCH
    
Number of semi-major axis used for the search of the range of SMA solving the resonance condition.
"""
const nbSMAResSEARCH = 1000

"""
    tabSMAResSEARCH

Table of semi-major axis, uniform in log-space, where we determine if the resonance condition can be solved in j-space.
ATTENTION, the table is sorted!
"""
const tabSMAResSEARCH = exp.(range(log(amin),length=nbSMAResSEARCH,log(amax))) 

##################################################

"""
    SMAResQ(np,ap,omega)
    
Checks if a resonance is possible for the semi-major axis ap with resonance number np.
Returns `true` if a given semi-major axis `ap` is such that there exists `jp` with `jlc(a)<jp<1.0` such that `np*nup(ap,jp) = omega`.
This function is very simple because jp->nup(ap,jp) is always decreasing.
"""
function SMAResQ(np::Int64,ap::Float64,omega::Float64)
    omegaleft = np*nup(ap,jlc(ap)) # Value of the resonant frequency on the very left
    omegaright = np*nup(ap,1.0) # Value of the resonant frequency on the very right
    return ( (omegaleft-omega)*(omegaright-omega) < 0.0) # If this product is negative, there is indeed a solution in the range jlc(ap)<j<1.0
end

"""
    getjRes(np,ap,omega)
    
Returns the value of `j` such that `np*nup(ap,jp) = omega`.
This function is very simple because jp->nup(ap,jp) is always decreasing.
This function will fail if the bisection condition is not satisfied: check with SMAResQ(np,ap,omega).
"""
function getjRes(np::Int64,ap::Float64,omega::Float64)
    return bisection(jp -> np*nup(ap,jp)-omega,jlc(ap),1.0) # We search for the root, by bisection, within the range jlc(ap)<j<1.0
end

"""
    tabSMARes!(n,np,a,j,tabSMARes_arg=tabSMARes)
    
Computes, if it exists, a given resonance line as an array of (ap,jp) locations in orbital space.
The values are written the the array tabSMARes_arg, which is set to be tabSMARes by default.

This function returns true if a line has been constructed and false otherwise.

# Remarks
- Part of the code could be refactored.
- The scan in SMAResSearch sometimes makes a lot of useless calculations.
"""
function tabSMARes!(n::Int64,np::Int64,a::Float64,j::Float64,tabSMARes_arg=tabSMARes)
    #####
    omega = n*nup(a,j) # Value of the resonant frequency to match
    #####
    # First, we scroll over tabSMAResSEARCH, to determine in that array
    # the smallest and largest SMA for which the resonance condition is solvable
    iSMAmin1, iSMAmax1 = nbSMAResSEARCH+1, 0 # Indices of the boundary of the SMA-range, taken outside the array range
    iSMAmin2, iSMAmax2 = nbSMAResSEARCH+1, 0 # Indices for the second resonance line (if any)

    ##### 
    take = -1 # Index to check which resonance line we are on
    nbLines = 0 # Number of resonance lines
    #####
    for iSMA = 1:nbSMAResSEARCH # Loop over the table of SMA used for the search of the SMA where the condition can be solved
        ap = tabSMAResSEARCH[iSMA] # Current value of ap for which we try to solve the resonance condition
        #####
        if(SMAResQ(np,ap,omega)) # The resonance condition can be solved for this value of ap
            #####
            if (take == -1) #enter first resonance line
                take = 0
                nbLines = 1
            elseif (take == 1) #enter second resonance line
                nbLines = 2
            end
            if (take == 0) #first resonance line
                if (iSMA < iSMAmin1) # We found a smaller value of ap for which the resonance condition can be solved
                    iSMAmin1 = iSMA # Updating iSMAmin  
                end
            elseif (take == 1) #second resonance line
                if (iSMA < iSMAmin2) # We found a smaller value of ap for which the resonance condition can be solved
                    iSMAmin2 = iSMA # Updating iSMAmin  
                end
            end
            #####
            if (take == 0) #first resonance line
                if (iSMA > iSMAmax1) # We found a larger value of ap for which the resonance condition can be solved
                    iSMAmax1 = iSMA # Updating iSMAmax
                end
            elseif (take == 1) #second resonance line
                if (iSMA > iSMAmax2) # We found a smaller value of ap for which the resonance condition can be solved
                    iSMAmax2 = iSMA # Updating iSMAmin  
                end
            end
        else
            if (take == 0) #end of the first (main) resonance line
                take = 1
            end
        end
    end
    #####
    # We have not been able to find SMA for which the resonance condition can be solved
    if (nbLines == 0)
        return 0
    end
    #####
    # In all the other cases, we have been able to find a SMA-region
    #####
    # First resonance line
    if (iSMAmin1 == 1) # The minimum SMA is exactly at the edge of the region
        tabSMARes_arg[1] = amin # Lower edge of the considered region
    #####
    else # The lower edge is not at the edge of the SMA domain
        # We must now find on which side of the (ap,jp)-triangle the resonant line ends, in the range given by tabSMAResSEARCH
        SMAlower, SMAupper = tabSMAResSEARCH[iSMAmin1-1], tabSMAResSEARCH[iSMAmin1] # Range in SMA over which we must search for a solution
        #####
        # First, we try to see if we can find a solution on the right, i.e. along jp=1.0
        omegalower, omegaupper = np*nup(SMAlower,1.0)-omega, np*nup(SMAupper,1.0)-omega # Value of the resonance frequency along jp=1.0, in the SMA range [SMAlower,SMAupper]
        #####
        if (omegalower*omegaupper < 0.0) # We found a bracketing interval along jp=1.0, i.e. the resonant location is on the right
            tabSMARes_arg[1] = bisection(ap -> np*nup(ap,1.0)-omega,SMAlower,SMAupper) # Finding the solution along jp=1.0, in the SMA range [SMAlower,SMAupper]
        else # Else, we must search for the resonant location along jp=jlc(ap), i.e. the resonant location is on the left
            tabSMARes_arg[1] = bisection(ap -> np*nup(ap,jlc(ap))-omega,SMAlower,SMAupper) # Finding the solution along jp=jlc(ap), in the SMA range [SMAlower,SMAupper]
        end
    end
    #####
    if (iSMAmax1 == nbSMAResSEARCH) # The maximum SMA is exactly at the edge of the region
        tabSMARes_arg[2] = amax # Upper edge of the considered region
    #####
    else # The upper edge is not at the edge of the SMA domain
        # We must now find on which side of the (ap,jp)-triangle the resonant line ends, in the range given by tabSMA ResSEARCH
        SMAlower, SMAupper = tabSMAResSEARCH[iSMAmax1], tabSMAResSEARCH[iSMAmax1+1] # Range in SMA over which we must search for a solution
        #####
        # First, we try to see if we can find a solution on the right, i.e. along jp=1.0
        omegalower, omegaupper = np*nup(SMAlower,1.0)-omega, np*nup(SMAupper,1.0)-omega # Value of the resonance frequency along jp=1.0, in the SMA range [SMAlower,SMAupper]
        #####
        if (omegalower*omegaupper < 0.0) # We found a bracketing interval along jp=1.0, i.e. the resonant location is on the right
            tabSMARes_arg[2] = bisection(ap -> np*nup(ap,1.0)-omega,SMAlower,SMAupper) # Finding the solution along jp=1.0, in the SMA range [SMAlower,SMAupper]
        else # Else, we must search for the resonant location along jp=jlc(ap), i.e. the resonant location is on the left
            tabSMARes_arg[2] = bisection(ap -> np*nup(ap,jlc(ap))-omega,SMAlower,SMAupper) # Finding the solution along jp=jlc(ap), in the SMA range [SMAlower,SMAupper]
        end
    end

    #####
    #####
    # Second resonance line (if any)
   
     if (nbLines == 2) # If there is asecond resonance line
    

        if (iSMAmin2 == 1) # The minimum SMA is exactly at the edge of the region
            tabSMARes_arg[3] = amin # Lower edge of the considered region
    #####
        else # The lower edge is not at the edge of the SMA domain
        # We must now find on which side of the (ap,jp)-triangle the resonant line ends, in the range given by tabSMAResSEARCH
            SMAlower, SMAupper = tabSMAResSEARCH[iSMAmin2-1], tabSMAResSEARCH[iSMAmin2] # Range in SMA over which we must search for a solution
        #####
        # First, we try to see if we can find a solution on the right, i.e. along jp=1.0
            omegalower, omegaupper = np*nup(SMAlower,1.0)-omega, np*nup(SMAupper,1.0)-omega # Value of the resonance frequency along jp=1.0, in the SMA range [SMAlower,SMAupper]

        #####
            if (omegalower*omegaupper < 0.0) # We found a bracketing interval along jp=1.0, i.e. the resonant location is on the right
                tabSMARes_arg[3] = bisection(ap -> np*nup(ap,1.0)-omega,SMAlower,SMAupper) # Finding the solution along jp=1.0, in the SMA range [SMAlower,SMAupper]
            else # Else, we must search for the resonant location along jp=jlc(ap), i.e. the resonant location is on the left
                tabSMARes_arg[3] = bisection(ap -> np*nup(ap,jlc(ap))-omega,SMAlower,SMAupper) # Finding the solution along jp=jlc(ap), in the SMA range [SMAlower,SMAupper]
            end
        end

    #####
        if (iSMAmax2 == nbSMAResSEARCH) # The maximum SMA is exactly at the edge of the region
            tabSMARes_arg[4] = amax # Upper edge of the considered region
    #####
        else # The upper edge is not at the edge of the SMA domain
        # We must now find on which side of the (ap,jp)-triangle the resonant line ends, in the range given by tabSMA ResSEARCH
            SMAlower, SMAupper = tabSMAResSEARCH[iSMAmax2], tabSMAResSEARCH[iSMAmax2+1] # Range in SMA over which we must search for a solution
        #####
        # First, we try to see if we can find a solution on the right, i.e. along jp=1.0
            omegalower, omegaupper = np*nup(SMAlower,1.0)-omega, np*nup(SMAupper,1.0)-omega # Value of the resonance frequency along jp=1.0, in the SMA range [SMAlower,SMAupper]
        #####
            if (omegalower*omegaupper < 0.0) # We found a bracketing interval along jp=1.0, i.e. the resonant location is on the right
                tabSMARes_arg[4] = bisection(ap -> np*nup(ap,1.0)-omega,SMAlower,SMAupper) # Finding the solution along jp=1.0, in the SMA range [SMAlower,SMAupper]
            else # Else, we must search for the resonant location along jp=jlc(ap), i.e. the resonant location is on the left
                tabSMARes_arg[4] = bisection(ap -> np*nup(ap,jlc(ap))-omega,SMAlower,SMAupper) # Finding the solution along jp=jlc(ap), in the SMA range [SMAlower,SMAupper]
            end
        end
    else # If there is only one resonance line, reset to 0 the SMA values of the second resonance line
        tabSMARes_arg[3] = 0.0
        tabSMARes_arg[4] = 0.0
    end
    #####
    return nbLines # We have been able to find the edge of the resonant SMA-domain
end
