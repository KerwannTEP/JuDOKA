"""
    Table_dKdJ
    
`IntTable` structure used in the computation of the j-gradient of the geometric factor Klnnp.
Use `Table_create!()` to initialize.
"""
const Table_dKdJ = Table_create!()

"""
    Table_dKdJp
    
`IntTable` structure used in the computation of the jp-gradient of the geometric factor Klnnp.
Use `Table_create!()` to initialize.
"""
const Table_dKdJp = Table_create!()

##################################################
# Functions used to compute the integrands of I1 and I2
##################################################
"""
    psi(j,cosf)
    
Radius over semi-major axis, at a given angular momentum and the cosine of a mean anomaly f.
"""
function psi(j::Float64,cosf::Float64) 
    # Reduced radius r/a
    return j^2/(1.0+sqrt(1.0-j^2)*cosf)
end

"""
    dpsidj(j,cosf)

j-derivative of the radius over semi-major axis, at a given angular momentum and the cosine of a mean anomaly f.
"""
function dpsidj(j::Float64,cosf::Float64) 
    # j-Derivative psi
    exc = sqrt(1.0-j^2)
    numerator = 2.0*j*(1+exc*cosf) + j^3/exc * cosf
    denumerator = (1.0 + exc*cosf)^2
    return numerator/denumerator
end

"""
    ddMdfdj(j,cosf)
    
j-derivative of the Jacobian, at a given angular momentum and the cosine of a true anomaly f.
"""
function ddMdfdj(j::Float64,cosf::Float64) 
    # j-Derivative of the Jacobian
    return (2.0*j*dpsidj(j,cosf)*psi(j,cosf) - psi(j,cosf)^2)/(j^2)
end

"""
    phi(l,a,j,cosf,rp)

Computes the min(r,rp)^(l)/max(r,rp)^(l+1) part of the integrand of the geometric factor K_nn'(a,j,a',j').

# Arguments
- `l   ::Int64`  : mode of the multipole expansion.
- `a   ::Float64`: semi-major axis of the first orbit.
- `j   ::Float64`: reduced angular momentum of the first orbit.
- `cosf::Float64`: cosine of the true anomaly of the first orbit.
- `rp  ::Float64`: radius of the second orbit.
"""
function phi(l::Int64,a::Float64,j::Float64,cosf::Float64,rp::Float64) 
    # Min/max function
    r = a*psi(j,cosf)
    return min(r,rp)^(l)/max(r,rp)^(l+1)
end

"""
    eta(l,a,j,cosf,rp)

Function used in the computation of dK/dj which returns l if r<r' and -(l+1) if r>r'.

# Arguments
- `l   ::Int64`  : mode of the multipole expansion.
- `a   ::Float64`: semi-major axis of the first orbit.
- `j   ::Float64`: reduced angular momentum of the first orbit.
- `cosf::Float64`: cosine of the true anomaly of the first orbit.
- `rp  ::Float64`: radius of the second orbit.
"""
function eta(l::Int64,a::Float64,j::Float64,cosf::Float64,rp::Float64) 
    # Eta function used in I2
    r = a*psi(j,cosf)
    if r < rp
        return l 
    else
        return -(l+1) 
    end
end

"""
    zeta(l,a,j,cosf,rp)

Integrand of the derivative of the integrand of the geometric factor K_nn'(a,j,a',j').

# Arguments
- `l   ::Int64`  : mode of the multipole expansion.
- `a   ::Float64`: semi-major axis of the first orbit.
- `j   ::Float64`: reduced angular momentum of the first orbit.
- `cosf::Float64`: cosine of the true anomaly of the first orbit.
- `rp  ::Float64`: radius of the second orbit.
"""
function zeta(l::Int64,a::Float64,j::Float64,cosf::Float64,rp::Float64) 
    # Function used in the integrand
    return ddMdfdj(j,cosf) + (psi(j,cosf))/j *eta(l,a,j,cosf,rp)*dpsidj(j,cosf)
end

##################################################
# Computing the geometric factor Klnnp and its derivatives.
##################################################
"""
    dKlnnpdj(l,n,np,a,jt,ap,jp,Table_dKdJ_arg=Table_dKdJ)

j-derivative (with respect to the first orbit) of the geometric factor K_nn'(a,j,a',j').

# Arguments
- `l   ::Int64`  : mode of the multipole expansion.
- `n   ::Int64`  : first resonance number of the orbit.
- `np  ::Int64`  : second resonance number of the orbit.
- `a   ::Float64`: semi-major axis of the first orbit.
- `jt  ::Float64`: reduced angular momentum of the first orbit.
- `ap  ::Float64`: semi-major axis of the second orbit.
- `jp  ::Float64`: reduced angular momentum of the second orbit.
- `Table_dKdJ_arg:: IntTable`: Structure which contains the tables needed for the computation of this integral.
"""
function dKlnnpdj(l::Int64,n::Int64,np::Int64,a::Float64,jt::Float64,ap::Float64,jp::Float64,Table_dKdJ_arg=Table_dKdJ)
    # Computing tabdeltadP and tabdeltadQ
    for j=1:nbK
        # Computing tabdeltadP
        x = 0.0
        for i=(Table_dKdJ_arg.tabw[j]+1):(Table_dKdJ_arg.tabw[j+1]) # ATTENTION, indices are shifted by one
            x += (Table_dKdJ_arg.tabh1[i]+l*Table_dKdJ_arg.tabh2[i])*(Table_dKdJ_arg.tabR[i]/Table_dKdJ_arg.tabRp[j])^(l)/(Table_dKdJ_arg.tabRp[j]) # Adding the contribution
        end

        Table_dKdJ_arg.tabdeltaP[j] = x # Updating tabdeltadP
        # Computing tabdeltadQ
        x = 0.0
        for i=(Table_dKdJ_arg.tabw[j+1]+1):(Table_dKdJ_arg.tabw[j+2]) # ATTENTION, indices are shifted by one
            x += (Table_dKdJ_arg.tabh1[i]-(l+1)*Table_dKdJ_arg.tabh2[i])*(Table_dKdJ_arg.tabRp[j]/Table_dKdJ_arg.tabR[i])^(l)/(Table_dKdJ_arg.tabR[i]) # Adding the contribution
        end
        Table_dKdJ_arg.tabdeltaQ[j] = x # Updating tabdeltadQ
    end
    # Computing the contributions from the sum dP
    x = Table_dKdJ_arg.tabdeltaP[1] # Initial value of tabdP[1]
    sumP = Table_dKdJ_arg.tabgp[1]*x # Initialising the contribution from dP with the value of tabdP[1]
    for j=1:(nbK-1) # Loop over the samples
        x = (Table_dKdJ_arg.tabRp[j]/Table_dKdJ_arg.tabRp[j+1])^(l+1)*x + Table_dKdJ_arg.tabdeltaP[j+1] # Computing the value of dP[j+1]
        sumP += Table_dKdJ_arg.tabgp[j+1]*x # Updating the total contribution from P with dP[j+1]
    end
    # Computing the contributions from the sum dQ
    x = Table_dKdJ_arg.tabdeltaQ[nbK] # Initial value of tabdQ[1]
    sumQ = Table_dKdJ_arg.tabgp[nbK]*x # Initialising the contribution from dQ with the value of tabdQ[nbK]
    for j=nbK:-1:2 # Loop over the samples. ATTENTION, this is a decreasing recurrence
        x = (Table_dKdJ_arg.tabRp[j-1]/Table_dKdJ_arg.tabRp[j])^(l)*x + Table_dKdJ_arg.tabdeltaQ[j-1] # Computing the value of dQ[j-1]
        sumQ += Table_dKdJ_arg.tabgp[j-1]*x # Updating the total contribution from dQ with dQ[j-1]
    end
    return (sumP+sumQ)/(nbK^(2)) # Computing the total contribution from both dP1 and dQ1
end

"""
    dKlnnpdjp(l,n,np,a,jt,ap,jp,Table_dKdJ_arg=Table_dKdJp)

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
- `Table_dKdJ_arg:: IntTable`: Structure which contains the tables needed for the computation of this integral.
"""
function dKlnnpdjp(l::Int64,n::Int64,np::Int64,a::Float64,jt::Float64,ap::Float64,jp::Float64,Table_dKdJ_arg=Table_dKdJp)
    return dKlnnpdj(l,np,n,ap,jp,a,jt,Table_dKdJ_arg)
end

"""
    A_dAnnpSQ(n,np,a,jt,ap,jp,Table_dKdJ_arg=Table_dKdJ,Table_dKdJ_arg=Table_dKdJp)

Computation of the coefficient AnnpSQ at (a,jt,a',j') and of its derivatives with respect to j (first orbit) and jp (second orbit).
Returns a triplet (AnnpSQ, dAnnpSQ/dj, dAnnpSQ/djp) in this order.

# Arguments
- `n   ::Int64`  : first resonance number of the orbit.
- `np  ::Int64`  : second resonance number of the orbit.
- `a   ::Float64`: semi-major axis of the first orbit.
- `jt  ::Float64`: reduced angular momentum of the first orbit.
- `ap  ::Float64`: semi-major axis of the second orbit.
- `jp  ::Float64`: reduced angular momentum of the second orbit.
- `Table_dKdJ_arg:: IntTable`: Structure which contains the tables needed for the computation of dK/dj and K.
- `Table_dKdJp_arg::IntTable`: Structure which contains the tables needed for the computation of dK/djp.

"""
function A_dAnnpSQ(n::Int64,np::Int64,a::Float64,j::Float64,ap::Float64,jp::Float64,Table_dKdJ_arg=Table_dKdJ,Table_dKdJp_arg=Table_dKdJp)
    # Getting rid of the cases where the coefficients is necessarily zero
    if (abs(n) > lmax || abs(np) > lmax || mod(n-np,2) != 0)
        return 0.0
    else # In that case the result is non-zero
        dAp = 0.0
        dA = 0.0 # Initialising the result
        A = 0.0
        
        Table_init!(Table_dKdJ_arg,n,np,a,j,ap,jp)
        Table_init!(Table_dKdJp_arg,np,n,ap,jp,a,j)
        
        for l = max(abs(n),abs(np)):2:lmax # Loop over the l that have non-vanishing contributions !! ATTENTION, the step of the loop is 2
            K = Klnnp(l,n,np,a,j,ap,jp,Table_dKdJ_arg)
            sgn = sign(K)
            dKp = dKlnnpdjp(l,n,np,a,j,ap,jp,Table_dKdJp_arg)
            dK = dKlnnpdj(l,n,np,a,j,ap,jp,Table_dKdJ_arg)
            dAp += (ylnSQ(l,n)*ylnSQ(l,np))/((2.0*l+1.0)^(3))*(abs(K)*sgn*dKp) # Adding the contribution from (l,n,np) !! ATTENTION, not to forget the squared
            dA += (ylnSQ(l,n)*ylnSQ(l,np))/((2.0*l+1.0)^(3))*(abs(K)*sgn*dK) # Adding the contribution from (l,n,np) !! ATTENTION, not to forget the squared
            A += (ylnSQ(l,n)*ylnSQ(l,np))/((2.0*l+1.0)^(3))*(K)^2 # Adding the contribution from (l,n,np) !! ATTENTION, not to forget the squared
        end
        dAp *= 32.0*pi^(2)
        dA *= 32.0*pi^(2) # Putting the prefactor
        A *= 16.0*pi^(2)
        return A, dA, dAp # Output
    end
end

##################################################
# Computing the diffusion coefficients from its multipole expansion.
##################################################

"""
    dintDJJnn(n,np,a,j,ap,jp,Table_dKdJ_arg=Table_dKdJ,Table_dKdJ_arg=Table_dKdJp)

Computation of the integrand of DRR_JJnn, contribution of the (n,n') resonance number to DRR, and its partial derivatives with respect to j (first orbit) and jp (second orbit).
Returns a triplet (dInt/dj, dInt/djp, Int) in this order, where Int is the integrand of DRR_JJnn.

# Arguments
- `n   ::Int64`  : first resonance number of the orbit.
- `np  ::Int64`  : second resonance number of the orbit.
- `a   ::Float64`: semi-major axis of the first orbit.
- `j   ::Float64`: reduced angular momentum of the first orbit.
- `ap  ::Float64`: semi-major axis of the second orbit.
- `jp  ::Float64`: reduced angular momentum of the second orbit.
- `Table_dKdJ_arg:: IntTable`: Structure which contains the tables needed for the computation of dK/dj and K.
- `Table_dKdJp_arg::IntTable`: Structure which contains the tables needed for the computation of dK/djp.

"""
function dintDJJnn(n::Int64,np::Int64,a::Float64,j::Float64,ap::Float64,jp::Float64,Table_dKdJ_arg=Table_dKdJ,Table_dKdJp_arg=Table_dKdJp)
    dnup = dnupdj(ap,jp)
    AdAdAp = A_dAnnpSQ(n,np,a,j,ap,jp,Table_dKdJ_arg,Table_dKdJp_arg)
    Ftot = FtotBath(ap,jp)
    # term dj'
    term1djp = NBathTot(ap) * (dfjBathdj(ap,jp)*abs(dnup) - fjBath(ap,jp)*d2nupdj2(ap,jp)*sign(dnup) )/(dnup)^2  * AdAdAp[1]
    term2djp = Ftot/abs(dnup) * AdAdAp[3]
    djp = term1djp + term2djp
    # term dj
    cst = Ftot/abs(dnup)
    dj = cst * AdAdAp[2]
    # term A
    termA = cst * AdAdAp[1]
    return dj, djp, termA
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
    DintDJJnnDj(n,np,a,j,ap,jp,Table_dKdJ_arg=Table_dKdJ,Table_dKdJ_arg=Table_dKdJp)

Computation of the integrand of DRR_JJnn, contribution of the (n,n') resonance number to DRR, and its total j-derivative.
Returns a couple (DInt/dj, Int) in this order, where Int is the integrand of DRR_JJnn.

# Arguments
- `n   ::Int64`  : first resonance number of the orbit.
- `np  ::Int64`  : second resonance number of the orbit.
- `a   ::Float64`: semi-major axis of the first orbit.
- `j   ::Float64`: reduced angular momentum of the first orbit.
- `ap  ::Float64`: semi-major axis of the second orbit.
- `jp  ::Float64`: reduced angular momentum of the second orbit.
- `Table_dKdJ_arg:: IntTable`: Structure which contains the tables needed for the computation of dK/dj and K.
- `Table_dKdJp_arg::IntTable`: Structure which contains the tables needed for the computation of dK/djp.

"""
function DintDJJnnDj(n::Int64,np::Int64,a::Float64,j::Float64,ap::Float64,jp::Float64,Table_dKdJ_arg=Table_dKdJ,Table_dKdJp_arg=Table_dKdJp)
    dDdj_dDdjp_termA = dintDJJnn(n,np,a,j,ap,jp,Table_dKdJ_arg,Table_dKdJp_arg)
    dj = djpdj(n,np,a,j,ap,jp)
    return dDdj_dDdjp_termA[1] + dj*dDdj_dDdjp_termA[2], dDdj_dDdjp_termA[3]
end

"""
    D_DRR_JJnnpDj(n,np,a,j,Table_dKdJ_arg=Table_dKdJ,Table_dKdJ_arg=Table_dKdJp)

Computation of DRR_JJnn, contribution of the (n,n') resonance number to DRR, and its total j-derivative.
Returns a couple (DRR_JJnn, D(DRR_JJnn)/Dj) in this order.

Using log-sampling in the integration for better interpolation.

# Arguments
- `n   ::Int64`  : first resonance number of the orbit.
- `np  ::Int64`  : second resonance number of the orbit.
- `a   ::Float64`: semi-major axis of the orbit.
- `j   ::Float64`: reduced angular momentum of the orbit.
- `Table_dKdJ_arg:: IntTable`: Structure which contains the tables needed for the computation of dK/dj and K.
- `Table_dKdJp_arg::IntTable`: Structure which contains the tables needed for the computation of dK/djp.

"""
function D_DRR_JJnnpDj(n::Int64,np::Int64,a::Float64,j::Float64,Table_dKdJ_arg=Table_dKdJ,Table_dKdJp_arg=Table_dKdJp,tabSMARes_arg=tabSMARes)
    omega = n*nup(a,j) # Value of the resonant frequency
    ResPoints = tabSMARes!(n,np,a,j,tabSMARes_arg) # Trying to compute the resonant line
    if (ResPoints == false)
        return 0.0, 0.0 # If the resonant line is empty, we return 0.0 for the diffusion coefficient
    end
    apmin, apmax = tabSMARes_arg[1], tabSMARes_arg[2] # Boundary of the SMA-region where the resonant line lies
    dalphap = log(apmax/apmin)/(nbResPoints) # Step distance in log-space
    dDnn = 0.0 # Initialisation of the result
    Dnn = 0.0
    termA, dinte = 0.0, 0.0
    for alphap = range(log(apmin),length=nbResPoints,log(apmax)-dalphap) # Uniform logarithmic spacing of the SMA !! ATTENTION, we do not go through the last element, as we will sum in the middle of the interval
        ap = exp(alphap + 0.5*dalphap) # Central position of the SMA considered
        jp = getjRes(np,ap,omega) # Finding the location of the resonant jp
        dinte, termA = DintDJJnnDj(n,np,a,j,ap,jp,Table_dKdJ_arg,Table_dKdJp_arg)
        Dnn += dalphap*ap*termA
        dDnn += dalphap*ap*dinte # Local contribution !! ATTENTION, not to forget the ap
    end
    Dnn *= (n^(2))/(abs(np)) # Adding the prefactor in n^2/|np|
    dDnn *= (n^(2))/(abs(np))
    return Dnn, dDnn#, apmin, apmax # Output
end

"""
    D_DRR_JJDj(a,j,Table_dKdJ_arg=Table_dKdJ,Table_dKdJ_arg=Table_dKdJp)

Computation of DRR_JJ the scalar resonant diffusion coefficient in JJ and its total j-derivative.
Returns a couple (DRR_JJ, D(DRR_JJ)/Dj) in this order.

# Arguments
- `a   ::Float64`: semi-major axis of the orbit.
- `j   ::Float64`: reduced angular momentum of the orbit.
- `Table_dKdJ_arg:: IntTable`: Structure which contains the tables needed for the computation of dK/dj and K.
- `Table_dKdJp_arg::IntTable`: Structure which contains the tables needed for the computation of dK/djp.

"""
function D_DRR_JJDj(a::Float64,j::Float64,Table_dKdJ_arg=Table_dKdJ,Table_dKdJp_arg=Table_dKdJp,tabSMARes_arg=tabSMARes)
    D = 0.0 # Initialisation of the result
    dD = 0.0
    for iPair=1:nbResPairs # Loop over the resonant pairs
        n, np = tabResPairs[1,iPair], tabResPairs[2,iPair] # Value of the considered resonance (n,np)
        Dnn_dDnn = D_DRR_JJnnpDj(n,np,a,j,Table_dKdJ_arg,Table_dKdJp_arg,tabSMARes_arg)
        D += Dnn_dDnn[1]
        dD += Dnn_dDnn[2]
    # Adding the contribution from the resonance pair (n,np)
    end
    D *= 4.0*pi*G^(2) # Adding the prefactor 4piG^2
    dD *= 4.0*pi*G^(2)
    return D, dD # Output
end


"""
    D_DRR_JJDJ(a,j,Table_dKdJ_arg=Table_dKdJ,Table_dKdJ_arg=Table_dKdJp)

Computation of DRR_JJ the scalar resonant diffusion coefficient in JJ and its total J-derivative.
Returns a couple (DRR_JJ, D(DRR_JJ)/DJ) in this order.

# Arguments
- `a   ::Float64`: semi-major axis of the orbit.
- `j   ::Float64`: reduced angular momentum of the orbit.
- `Table_dKdJ_arg:: IntTable`: Structure which contains the tables needed for the computation of dK/dj and K.
- `Table_dKdJp_arg::IntTable`: Structure which contains the tables needed for the computation of dK/djp.

"""
function D_DRR_JJDJ(a::Float64,j::Float64,Table_dKdJ_arg=Table_dKdJ,Table_dKdJp_arg=Table_dKdJp,tabSMARes_arg=tabSMARes)
    D_dD = D_DRR_JJDj(a,j,Table_dKdJ_arg,Table_dKdJp_arg,tabSMARes_arg)
    return D_dD[1], D_dD[2]/Jc(a)
end


"""
    DRR_J_JJ(a,j,Table_dKdJ_arg=Table_dKdJ,Table_dKdJ_arg=Table_dKdJp)

Computation of DRR_JJ the scalar resonant diffusion coefficient in JJ-coordinates and of DRR_J the scalar resonant diffusion coefficient in J-coordinate.
Returns a couple (DRR_J, DRR_JJ) in this order.

# Arguments
- `a   ::Float64`: semi-major axis of the orbit.
- `j   ::Float64`: reduced angular momentum of the orbit.
- `Table_dKdJ_arg:: IntTable`: Structure which contains the tables needed for the computation of dK/dj and K.
- `Table_dKdJp_arg::IntTable`: Structure which contains the tables needed for the computation of dK/djp.

"""
function DRR_J_JJ(a::Float64,j::Float64,Table_dKdJ_arg=Table_dKdJ,Table_dKdJp_arg=Table_dKdJp,tabSMARes_arg=tabSMARes)
    DrrJJ_dDrrJJ = D_DRR_JJDJ(a,j,Table_dKdJ_arg,Table_dKdJp_arg,tabSMARes_arg)
    J = j*Jc(a)
    return DrrJJ_dDrrJJ[1]/(2.0*J) + DrrJJ_dDrrJJ[2]/2, DrrJJ_dDrrJJ[1] # Returns DRR_J, DRR_JJ
end

"""
    DRR_J(a,j,Table_dKdJ_arg=Table_dKdJ,Table_dKdJ_arg=Table_dKdJp)

Computation of DRR_J the scalar resonant diffusion coefficient in J-coordinate.

# Arguments
- `a   ::Float64`: semi-major axis of the orbit.
- `j   ::Float64`: reduced angular momentum of the orbit.
- `Table_dKdJ_arg:: IntTable`: Structure which contains the tables needed for the computation of dK/dj and K.
- `Table_dKdJp_arg::IntTable`: Structure which contains the tables needed for the computation of dK/djp.

"""
function DRR_J(a::Float64,j::Float64,Table_dKdJ_arg=Table_dKdJ,Table_dKdJp_arg=Table_dKdJp,tabSMARes_arg=tabSMARes)
    return DRR_J_JJ(a,j,Table_dKdJ_arg,Table_dKdJp_arg,tabSMARes_arg)[1] # Returns DRR_J
end

##################################################
# Recovering the diffusion coefficients in j and jj.
##################################################

"""
    DRR_j_jj(a,j,Table_dKdJ_arg=Table_dKdJ,Table_dKdJ_arg=Table_dKdJp)

Computation of DRR_jj the (normalized) scalar resonant diffusion coefficient in jj-coordinates and of DRR_j the scalar resonant diffusion coefficient in j-coordinate.
Returns a couple (DRR_j, DRR_jj) in this order.

# Arguments
- `a   ::Float64`: semi-major axis of the orbit.
- `j   ::Float64`: reduced angular momentum of the orbit.
- `Table_dKdJ_arg:: IntTable`: Structure which contains the tables needed for the computation of dK/dj and K.
- `Table_dKdJp_arg::IntTable`: Structure which contains the tables needed for the computation of dK/djp.

"""
function DRR_j_jj(a::Float64,j::Float64,Table_dKdJ_arg=Table_dKdJ,Table_dKdJp_arg=Table_dKdJp,tabSMARes_arg=tabSMARes)
    DrrJ_DrrJJ = DRR_J_JJ(a,j,Table_dKdJ_arg,Table_dKdJp_arg,tabSMARes_arg)
    return DrrJ_DrrJJ[1]/Jc(a), DrrJ_DrrJJ[2]/(Jc(a)^(2))# Returns DRR_j, DRR_jj
end

"""
    DRR_j(a,j,Table_dKdJ_arg=Table_dKdJ,Table_dKdJ_arg=Table_dKdJp)

Computation of DRR_j the (normalized) scalar resonant diffusion coefficient in j-coordinate.

# Arguments
- `a   ::Float64`: semi-major axis of the orbit.
- `j   ::Float64`: reduced angular momentum of the orbit.
- `Table_dKdJ_arg:: IntTable`: Structure which contains the tables needed for the computation of dK/dj and K.
- `Table_dKdJp_arg::IntTable`: Structure which contains the tables needed for the computation of dK/djp.

"""
function DRR_j(a::Float64,j::Float64,Table_dKdJ_arg=Table_dKdJ,Table_dKdJp_arg=Table_dKdJp,tabSMARes_arg=tabSMARes)    
    return DRR_J(a,j,Table_dKdJ_arg,Table_dKdJp_arg,tabSMARes_arg)[1]/(Jc(a))# Returns DRR_j
end