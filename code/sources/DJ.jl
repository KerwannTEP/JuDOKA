"""
    IntTable_dKdJ
    
`IntTable` structure used in the computation of the j-gradient of the geometric factor Klnnp.
Use `Table_create!()` to initialize.
"""
const IntTable_dKdJ_serial = IntTable_create!()

"""
    IntTable_dKdJp
    
`IntTable` structure used in the computation of the jp-gradient of the geometric factor Klnnp.
Use `Table_create!()` to initialize.
"""
const IntTable_dKdJp_serial = IntTable_create!()

##################################################
# Functions used to compute the integrands of I1 and I2
##################################################
"""
    reduced_rad(j,cosf)
    
Radius over semi-major axis, at a given angular momentum and the cosine of a mean anomaly f.
"""
function reduced_rad(j::Float64,cosf::Float64) 
    # Reduced radius r/a
    return j^2/(1.0+sqrt(1.0-j^2)*cosf)
end

"""
    Dreduced_rad(j,cosf)

j-derivative of the radius over semi-major axis, at a given angular momentum and the cosine of a mean anomaly f.
"""
function Dreduced_rad(j::Float64,cosf::Float64) 
    # j-Derivative psi
    ecc = sqrt(1.0-j^2)
    numerator   = 2.0*j*(1+ecc*cosf) + j^3/ecc * cosf
    denumerator = (1.0 + ecc*cosf)^2
    return numerator/denumerator
end

"""
    ddMdfdj(j,cosf)
    
j-derivative of the Jacobian, at a given angular momentum and the cosine of a true anomaly f.
"""
function ddMdfdj(j::Float64,cosf::Float64) 
    # j-Derivative of the Jacobian
    return (2.0*j*Dreduced_rad(j,cosf)*reduced_rad(j,cosf) - reduced_rad(j,cosf)^2)/(j^2)
end

##################################################
# Computing the geometric factor Klnnp and its derivatives.
##################################################
"""
    dKlnnpdj(l,n,np,a,jt,ap,jp,IntTable_dKdJ=IntTable_dKdJ_serial)

j-derivative (with respect to the first orbit) of the geometric factor K_nn'(a,j,a',j').

# Arguments
- `l   ::Int64`  : mode of the multipole expansion.
- `n   ::Int64`  : first resonance number of the orbit.
- `np  ::Int64`  : second resonance number of the orbit.
- `a   ::Float64`: semi-major axis of the first orbit.
- `jt  ::Float64`: reduced angular momentum of the first orbit.
- `ap  ::Float64`: semi-major axis of the second orbit.
- `jp  ::Float64`: reduced angular momentum of the second orbit.
- `IntTable_dKdJ :: IntTable`: Structure which contains the tables needed for the computation of this integral.
"""
function dKlnnpdj(l::Int64,n::Int64,np::Int64,a::Float64,jt::Float64,ap::Float64,jp::Float64,
                  IntTable_dKdJ::IntTable=IntTable_dKdJ_serial)

    tabw  = IntTable_dKdJ.tabw
    tabgp = IntTable_dKdJ.tabgp
    tabh1 = IntTable_dKdJ.tabh1
    tabh2 = IntTable_dKdJ.tabh2
    tabR  = IntTable_dKdJ.tabR
    tabRp = IntTable_dKdJ.tabRp
    tabdeltaP = IntTable_dKdJ.tabdeltaP
    tabdeltaQ = IntTable_dKdJ.tabdeltaQ

    # Computing tabdeltadP and tabdeltadQ
    for j=1:nbK
        # Computing tabdeltadP
        x = 0.0
        for i=(tabw[j]+1):(tabw[j+1]) # ATTENTION, indices are shifted by one
            x += (tabh1[i]+l*tabh2[i])*(tabR[i]/tabRp[j])^(l)/(tabRp[j]) # Adding the contribution
        end

        tabdeltaP[j] = x # Updating tabdeltadP
        # Computing tabdeltadQ
        x = 0.0
        for i=(tabw[j+1]+1):(tabw[j+2]) # ATTENTION, indices are shifted by one
            x += (tabh1[i]-(l+1)*tabh2[i])*(tabRp[j]/tabR[i])^(l)/(tabR[i]) # Adding the contribution
        end
        tabdeltaQ[j] = x # Updating tabdeltadQ
    end
    # Computing the contributions from the sum dP
    x = tabdeltaP[1] # Initial value of tabdP[1]
    sumP = tabgp[1]*x # Initialising the contribution from dP with the value of tabdP[1]
    for j=1:(nbK-1) # Loop over the samples
        x = (tabRp[j]/tabRp[j+1])^(l+1)*x + tabdeltaP[j+1] # Computing the value of dP[j+1]
        sumP += tabgp[j+1]*x # Updating the total contribution from P with dP[j+1]
    end
    # Computing the contributions from the sum dQ
    x = tabdeltaQ[nbK] # Initial value of tabdQ[1]
    sumQ = tabgp[nbK]*x # Initialising the contribution from dQ with the value of tabdQ[nbK]
    for j=nbK:-1:2 # Loop over the samples. ATTENTION, this is a decreasing recurrence
        x = (tabRp[j-1]/tabRp[j])^(l)*x + tabdeltaQ[j-1] # Computing the value of dQ[j-1]
        sumQ += tabgp[j-1]*x # Updating the total contribution from dQ with dQ[j-1]
    end
    return (sumP+sumQ)/(nbK^(2)) # Computing the total contribution from both dP1 and dQ1
end

"""
    dKlnnpdjp(l,n,np,a,jt,ap,jp,IntTable_dKdJ=IntTable_dKdJp_serial)

jp-derivative (with respect to the second orbit) of the geometric factor K_nn'(a,j,a',j').
Uses the relation dK/dj'(n,n',a,j,a',j') = dK/dj(n',n,a',j',a,j).

# Arguments
- `l   ::Int64`  : mode of the multipole expansion.
- `n   ::Int64`  : first resonance number of the orbit.
- `np  ::Int64`  : second resonance number of the orbit.
- `a   ::Float64`: semi-major axis of the first orbit.
- `jt  ::Float64`: reduced angular momentum of the first orbit.
- `ap  ::Float64`: semi-major axis of the second orbit.
- `jp  ::Float64`: reduced angular momentum of the second orbit.
- `IntTable_dKdJ:: IntTable`: Structure which contains the tables needed for the computation of this integral.
"""
function dKlnnpdjp(l::Int64,n::Int64,np::Int64,a::Float64,jt::Float64,ap::Float64,jp::Float64,
                   IntTable_dKdJ::IntTable=IntTable_dKdJp_serial)

    return dKlnnpdj(l,np,n,ap,jp,a,jt,IntTable_dKdJ)
end

"""
    A_dAnnpSQ(n,np,a,jt,ap,jp,IntTable_dKdJ=IntTable_dKdJ_serial,IntTable_dKdJ=IntTable_dKdJp_serial)

Computation of the coefficient AnnpSQ at (a,jt,a',j') and of its derivatives with respect to j (first orbit) and jp (second orbit).
Returns a triplet (AnnpSQ, dAnnpSQ/dj, dAnnpSQ/djp) in this order.

# Arguments
- `n   ::Int64`  : first resonance number of the orbit.
- `np  ::Int64`  : second resonance number of the orbit.
- `a   ::Float64`: semi-major axis of the first orbit.
- `jt  ::Float64`: reduced angular momentum of the first orbit.
- `ap  ::Float64`: semi-major axis of the second orbit.
- `jp  ::Float64`: reduced angular momentum of the second orbit.
- `IntTable_dKdJ:: IntTable`: Structure which contains the tables needed for the computation of dK/dj and K.
- `IntTable_dKdJp::IntTable`: Structure which contains the tables needed for the computation of dK/djp.

"""
function A_dAnnpSQ(n::Int64,np::Int64,a::Float64,j::Float64,ap::Float64,jp::Float64,
                   IntTable_dKdJ::IntTable=IntTable_dKdJ_serial,IntTable_dKdJp::IntTable=IntTable_dKdJp_serial)

    # Getting rid of the cases where the coefficients is necessarily zero
    if (abs(n) > lmax || abs(np) > lmax || mod(n-np,2) != 0)
        return 0.0
    else # In that case the result is non-zero

        # Initialising the result
        A     = 0.0
        dAdj  = 0.0
        dAdjp = 0.0
        
        IntTable_init!(IntTable_dKdJ ,n ,np,a ,j ,ap,jp)
        IntTable_init!(IntTable_dKdJp,np,n ,ap,jp,a ,j)
        
        for l = max(abs(n),abs(np)):2:lmax # Loop over the l that have non-vanishing contributions !! ATTENTION, the step of the loop is 2
            K   =     Klnnp(l,n,np,a,j,ap,jp,IntTable_dKdJ)
            dK  =  dKlnnpdj(l,n,np,a,j,ap,jp,IntTable_dKdJ)
            dKp = dKlnnpdjp(l,n,np,a,j,ap,jp,IntTable_dKdJp)
            sgn = sign(K)
            
            A     += (ylnSQ(l,n)*ylnSQ(l,np))/((2.0*l+1.0)^(3))*(K)^2 
            dAdj  += (ylnSQ(l,n)*ylnSQ(l,np))/((2.0*l+1.0)^(3))*(abs(K)*sgn*dK) 
            dAdjp += (ylnSQ(l,n)*ylnSQ(l,np))/((2.0*l+1.0)^(3))*(abs(K)*sgn*dKp) 
        end

        # Putting the prefactor
        A     *= 16.0*pi^(2)
        dAdj  *= 32.0*pi^(2)
        dAdjp *= 32.0*pi^(2)
        
        return A, dAdj, dAdjp # Output
    end
end

##################################################
# Computing the diffusion coefficients from its multipole expansion.
##################################################

"""
    dintDJJnn(n,np,a,j,ap,jp,IntTable_dKdJ=IntTable_dKdJ_serial,IntTable_dKdJ=IntTable_dKdJp_serial)

Computation of the integrand of DRR_JJnn, contribution of the (n,n') resonance number to DRR, and its partial derivatives with respect to j (first orbit) and jp (second orbit).
Returns a triplet (Int,dInt/dj, dInt/djp) in this order, where Int is the integrand of DRR_JJnn.

# Arguments
- `n   ::Int64`  : first resonance number of the orbit.
- `np  ::Int64`  : second resonance number of the orbit.
- `a   ::Float64`: semi-major axis of the first orbit.
- `j   ::Float64`: reduced angular momentum of the first orbit.
- `ap  ::Float64`: semi-major axis of the second orbit.
- `jp  ::Float64`: reduced angular momentum of the second orbit.
- `IntTable_dKdJ:: IntTable`: Structure which contains the tables needed for the computation of dK/dj and K.
- `IntTable_dKdJp::IntTable`: Structure which contains the tables needed for the computation of dK/djp.

"""
function dintDJJnn(n::Int64,np::Int64,a::Float64,j::Float64,ap::Float64,jp::Float64,
                   IntTable_dKdJ::IntTable=IntTable_dKdJ_serial,IntTable_dKdJp::IntTable=IntTable_dKdJp_serial)

    dnup = dnupdj(ap,jp)
    A, dAdj, dAdjp = A_dAnnpSQ(n,np,a,j,ap,jp,IntTable_dKdJ,IntTable_dKdJp)
    Ftot = FtotBath(ap,jp)
    cst = Ftot/abs(dnup)

    # Integrand of DRR_JJnn
    Integrand = cst * A

    # j-gradient of the integrand of DRR_JJnn
    dIntdj = cst * dAdj

    # j'-gradient of the integrand of DRR_JJnn
    term1djp = NBathTot(ap) * (dfjBathdj(ap,jp)*abs(dnup) - fjBath(ap,jp)*d2nupdj2(ap,jp)*sign(dnup) )/(dnup)^2  * A
    term2djp = Ftot/abs(dnup) * dAdjp
    dIntdjp = term1djp + term2djp
   
    return Integrand, dIntdj, dIntdjp
end

"""
    djpdj(n,np,a,j,ap,jp)

Total j-derivative of jp(j), reduced angular momentum of the second orbit.
Uses implicit differentiation of the resonance delta-condition n*nu_P(a,j) = np*nu_P(a',j') where nu_P is the precession frequency.

# Arguments
- `n   ::Int64`  : first resonance number of the orbit.
- `np  ::Int64`  : second resonance number of the orbit.
- `a   ::Float64`: semi-major axis of the first orbit.
- `j   ::Float64`: reduced angular momentum of the first orbit.
- `ap  ::Float64`: semi-major axis of the second orbit.
- `jp  ::Float64`: reduced angular momentum of the second orbit.

"""
function djpdj(n::Int64,np::Int64,a::Float64,j::Float64,ap::Float64,jp::Float64)
    return (n*dnupdj(a,j))/(np*dnupdj(ap,jp))
end

"""
    DintDJJnnDj(n,np,a,j,ap,jp,IntTable_dKdJ=IntTable_dKdJ_serial,IntTable_dKdJ=IntTable_dKdJp_serial)

Computation of the integrand of DRR_JJnn, contribution of the (n,n') resonance number to DRR, and its total j-derivative.
Returns a couple (Int, DInt/Dj) in this order, where Int is the integrand of DRR_JJnn.

# Arguments
- `n   ::Int64`  : first resonance number of the orbit.
- `np  ::Int64`  : second resonance number of the orbit.
- `a   ::Float64`: semi-major axis of the first orbit.
- `j   ::Float64`: reduced angular momentum of the first orbit.
- `ap  ::Float64`: semi-major axis of the second orbit.
- `jp  ::Float64`: reduced angular momentum of the second orbit.
- `IntTable_dKdJ:: IntTable`: Structure which contains the tables needed for the computation of dK/dj and K.
- `IntTable_dKdJp::IntTable`: Structure which contains the tables needed for the computation of dK/djp.

"""
function DintDJJnnDj(n::Int64,np::Int64,a::Float64,j::Float64,ap::Float64,jp::Float64,
                     IntTable_dKdJ::IntTable=IntTable_dKdJ_serial,IntTable_dKdJp::IntTable=IntTable_dKdJp_serial)

    Integrand, dIntdj, dIntdjp = dintDJJnn(n,np,a,j,ap,jp,IntTable_dKdJ,IntTable_dKdJp)
    dj = djpdj(n,np,a,j,ap,jp)
    return Integrand, dIntdj + dj*dIntdjp
end

"""
    D_DRR_JJnnpDj(n,np,a,j,IntTable_dKdJ=IntTable_dKdJ_serial,IntTable_dKdJ=IntTable_dKdJp_serial,tabSMARes=tabSMARes_serial)

Computation of DRR_JJnn, contribution of the (n,n') resonance number to DRR, and its total j-derivative.
Returns a couple (DRR_JJnn, D(DRR_JJnn)/Dj) in this order.

Using log-sampling in the integration for better interpolation.

# Arguments
- `n   ::Int64`  : first resonance number of the orbit.
- `np  ::Int64`  : second resonance number of the orbit.
- `a   ::Float64`: semi-major axis of the orbit.
- `j   ::Float64`: reduced angular momentum of the orbit.
- `IntTable_dKdJ:: IntTable`: Structure which contains the tables needed for the computation of dK/dj and K.
- `IntTable_dKdJp::IntTable`: Structure which contains the tables needed for the computation of dK/djp.

"""
function D_DRR_JJnnpDj(n::Int64,np::Int64,a::Float64,j::Float64,IntTable_dKdJ::IntTable=IntTable_dKdJ_serial,
                       IntTable_dKdJp::IntTable=IntTable_dKdJp_serial,tabSMARes::Array{Float64,2}=tabSMARes_serial)

    omega = n*nup(a,j) # Value of the resonant frequency
    ResPoints = tabSMARes!(n,np,a,j,tabSMARes) # Trying to compute the resonant line
    if (ResPoints == 0)
        return 0.0, 0.0 # If the resonant line is empty, we return 0.0 for the diffusion coefficient
    end

    # First resonance line
    apmin1, apmax1 = tabSMARes[1,1], tabSMARes[1,2] # Boundary of the SMA-region where the resonant line lies
    dalphap1 = log(apmax1/apmin1)/(nbResPoints) # Step distance in log-space

    if (ResPoints == 2) # If second resonance line
        apmin2, apmax2 = tabSMARes[2,1], tabSMARes[2,2] # Boundary of the SMA-region where the resonant line lies
        dalphap2 = log(apmax2/apmin2)/(nbResPoints) # Step distance in log-space
    end

    dDnn = 0.0 # Initialisation of the result
    Dnn  = 0.0

    for alphap1 = range(log(apmin1),length=nbResPoints,log(apmax1)-dalphap1) # Uniform logarithmic spacing of the SMA !! ATTENTION, we do not go through the last element, as we will sum in the middle of the interval
        ap = exp(alphap1 + 0.5*dalphap1) # Central position of the SMA considered
        jp = getjRes(np,ap,omega) # Finding the location of the resonant jp
        Integrand, DIntDj = DintDJJnnDj(n,np,a,j,ap,jp,IntTable_dKdJ,IntTable_dKdJp)
        Dnn  += dalphap1*ap*Integrand
        dDnn += dalphap1*ap*DIntDj # Local contribution !! ATTENTION, not to forget the ap
    end

    if (ResPoints == 2) # If two resonance lines
        for alphap2 = range(log(apmin2),length=nbResPoints,log(apmax2)-dalphap2) # Uniform logarithmic spacing of the SMA !! ATTENTION, we do not go through the last element, as we will sum in the middle of the interval
            ap = exp(alphap2 + 0.5*dalphap2) # Central position of the SMA considered
            jp = getjRes(np,ap,omega) # Finding the location of the resonant jp
            Integrand, DIntDj = DintDJJnnDj(n,np,a,j,ap,jp,IntTable_dKdJ,IntTable_dKdJp)
            Dnn  += dalphap2*ap*Integrand
            dDnn += dalphap2*ap*DIntDj # Local contribution !! ATTENTION, not to forget the ap
        end
    end

    Dnn  *= (n^(2))/(abs(np)) # Adding the prefactor in n^2/|np|
    dDnn *= (n^(2))/(abs(np))
    return Dnn, dDnn#, apmin, apmax # Output
end

"""
    D_DRR_JJDj(a,j,IntTable_dKdJ=IntTable_dKdJ_serial,IntTable_dKdJ=IntTable_dKdJp_serial,tabSMARes=tabSMARes_serial)

Computation of DRR_JJ the scalar resonant diffusion coefficient in JJ and its total j-derivative.
Returns a couple (DRR_JJ, D(DRR_JJ)/Dj) in this order.

# Arguments
- `a   ::Float64`: semi-major axis of the orbit.
- `j   ::Float64`: reduced angular momentum of the orbit.
- `IntTable_dKdJ:: IntTable`: Structure which contains the tables needed for the computation of dK/dj and K.
- `IntTable_dKdJp::IntTable`: Structure which contains the tables needed for the computation of dK/djp.

"""
function D_DRR_JJDj(a::Float64,j::Float64,IntTable_dKdJ::IntTable=IntTable_dKdJ_serial,
                    IntTable_dKdJp::IntTable=IntTable_dKdJp_serial,tabSMARes::Array{Float64,2}=tabSMARes_serial)

    D  = 0.0 # Initialisation of the result
    dD = 0.0
    for iPair=1:nbResPairs # Loop over the resonant pairs
        n, np = tabResPairs[1,iPair], tabResPairs[2,iPair] # Value of the considered resonance (n,np)
        Dnn, dDnn = D_DRR_JJnnpDj(n,np,a,j,IntTable_dKdJ,IntTable_dKdJp,tabSMARes)
        D  += Dnn
        dD += dDnn
    # Adding the contribution from the resonance pair (n,np)
    end
    D  *= 4.0*pi*G^(2) # Adding the prefactor 4piG^2
    dD *= 4.0*pi*G^(2)
    return D, dD # Output
end


"""
    D_DRR_JJDJ(a,j,IntTable_dKdJ=IntTable_dKdJ_serial,IntTable_dKdJ=IntTable_dKdJp_serial,tabSMARes=tabSMARes_serial)

Computation of DRR_JJ the scalar resonant diffusion coefficient in JJ and its total J-derivative.
Returns a couple (DRR_JJ, D(DRR_JJ)/DJ) in this order.

# Arguments
- `a   ::Float64`: semi-major axis of the orbit.
- `j   ::Float64`: reduced angular momentum of the orbit.
- `IntTable_dKdJ:: IntTable`: Structure which contains the tables needed for the computation of dK/dj and K.
- `IntTable_dKdJp::IntTable`: Structure which contains the tables needed for the computation of dK/djp.

"""
function D_DRR_JJDJ(a::Float64,j::Float64,IntTable_dKdJ::IntTable=IntTable_dKdJ_serial,
                    IntTable_dKdJp::IntTable=IntTable_dKdJp_serial,tabSMARes::Array{Float64,2}=tabSMARes_serial)

    D, dD = D_DRR_JJDj(a,j,IntTable_dKdJ,IntTable_dKdJp,tabSMARes)
    return D, dD/Jc(a)
end


"""
    DRR_J_JJ(a,j,IntTable_dKdJ=IntTable_dKdJ_serial,IntTable_dKdJp=IntTable_dKdJp_serial,tabSMARes=tabSMARes_serial)

Computation of DRR_JJ the scalar resonant diffusion coefficient in JJ-coordinates and of DRR_J the scalar resonant diffusion coefficient in J-coordinate.
Returns a couple (DRR_J, DRR_JJ) in this order.

# Arguments
- `a   ::Float64`: semi-major axis of the orbit.
- `j   ::Float64`: reduced angular momentum of the orbit.
- `IntTable_dKdJ:: IntTable`: Structure which contains the tables needed for the computation of dK/dj and K.
- `IntTable_dKdJp::IntTable`: Structure which contains the tables needed for the computation of dK/djp.

"""
function DRR_J_JJ(a::Float64,j::Float64,IntTable_dKdJ::IntTable=IntTable_dKdJ_serial,
                  IntTable_dKdJp::IntTable=IntTable_dKdJp_serial,tabSMARes::Array{Float64,2}=tabSMARes_serial)

    DrrJJ, dDrrJJ = D_DRR_JJDJ(a,j,IntTable_dKdJ,IntTable_dKdJp,tabSMARes)
    J = j*Jc(a)
    return DrrJJ/(2.0*J) + dDrrJJ/2, DrrJJ # Returns DRR_J, DRR_JJ
end

"""
    DRR_J(a,j,IntTable_dKdJ=IntTable_dKdJ_serial,IntTable_dKdJ=IntTable_dKdJp_serial,tabSMARes=tabSMARes_serial)

Computation of DRR_J the scalar resonant diffusion coefficient in J-coordinate.

# Arguments
- `a   ::Float64`: semi-major axis of the orbit.
- `j   ::Float64`: reduced angular momentum of the orbit.
- `IntTable_dKdJ:: IntTable`: Structure which contains the tables needed for the computation of dK/dj and K.
- `IntTable_dKdJp::IntTable`: Structure which contains the tables needed for the computation of dK/djp.

"""
function DRR_J(a::Float64,j::Float64,IntTable_dKdJ::IntTable=IntTable_dKdJ_serial,
               IntTable_dKdJp::IntTable=IntTable_dKdJp_serial,tabSMARes::Array{Float64,2}=tabSMARes_serial)

    return DRR_J_JJ(a,j,IntTable_dKdJ,IntTable_dKdJp,tabSMARes)[1] # Returns DRR_J
end

##################################################
# Recovering the diffusion coefficients in j and jj.
##################################################

"""
    DRR_j_jj(a,j,IntTable_dKdJ=IntTable_dKdJ_serial,IntTable_dKdJ=IntTable_dKdJp_serial,tabSMARes=tabSMARes_serial)

Computation of DRR_jj the (normalized) scalar resonant diffusion coefficient in jj-coordinates and of DRR_j the scalar resonant diffusion coefficient in j-coordinate.
Returns a couple (DRR_j, DRR_jj) in this order.

# Arguments
- `a   ::Float64`: semi-major axis of the orbit.
- `j   ::Float64`: reduced angular momentum of the orbit.
- `IntTable_dKdJ:: IntTable`: Structure which contains the tables needed for the computation of dK/dj and K.
- `IntTable_dKdJp::IntTable`: Structure which contains the tables needed for the computation of dK/djp.

"""
function DRR_j_jj(a::Float64,j::Float64,IntTable_dKdJ::IntTable=IntTable_dKdJ_serial,
                  IntTable_dKdJp::IntTable=IntTable_dKdJp_serial,tabSMARes::Array{Float64,2}=tabSMARes_serial)

    DrrJ, DrrJJ = DRR_J_JJ(a,j,IntTable_dKdJ,IntTable_dKdJp,tabSMARes)
    return DrrJ/Jc(a), DrrJJ/(Jc(a)^(2))# Returns DRR_j, DRR_jj
end

"""
    DRR_j(a,j,IntTable_dKdJ=IntTable_dKdJ_serial,IntTable_dKdJ=IntTable_dKdJp_serial,tabSMARes=tabSMARes_serial)

Computation of DRR_j the (normalized) scalar resonant diffusion coefficient in j-coordinate.

# Arguments
- `a   ::Float64`: semi-major axis of the orbit.
- `j   ::Float64`: reduced angular momentum of the orbit.
- `IntTable_dKdJ:: IntTable`: Structure which contains the tables needed for the computation of dK/dj and K.
- `IntTable_dKdJp::IntTable`: Structure which contains the tables needed for the computation of dK/djp.

"""
function DRR_j(a::Float64,j::Float64,IntTable_dKdJ::IntTable=IntTable_dKdJ_serial,
               IntTable_dKdJp::IntTable=IntTable_dKdJp_serial,tabSMARes::Array{Float64,2}=tabSMARes_serial)  
  
    return DRR_J(a,j,IntTable_dKdJ,IntTable_dKdJp,tabSMARes)[1]/(Jc(a))# Returns DRR_j
end