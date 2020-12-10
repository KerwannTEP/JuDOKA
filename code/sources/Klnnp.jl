"""
    tabf_Klnnp
    
Table of the true anomalies, f, used for the discretisation of the integrals in the geometric factor Klnnp.
Set as `SVector{nbK,Float64}([(pi/nbK)*(k-0.5) for k=1:nbK])` by default.
Using `StaticArrays` for better performances.
"""
const tabf_Klnnp = SVector{nbK,Float64}([(pi/nbK)*(k-0.5) for k=1:nbK])

"""
    tabcosf_Klnnp
    
Table of cos(f) used to compute the radius in the geometric factor Klnnp, where f are the true anomalies.
Set as `SVector{nbK,Float64}([cos(tabf_Klnnp[k]) for k=1:nbK])` by default.
Using `StaticArrays` for better performances.
"""
const tabcosf_Klnnp = SVector{nbK,Float64}([cos(tabf_Klnnp[k]) for k=1:nbK]) 

"""
    tabcosnf_Klnnp
    
Table of cos(nf) used to compute the integrand in Klnnp the geometric factor Klnnp, where f are the true anomalies.
Set as `SMatrix{lmax,nbK,Float64}([cos(n*tabf_Klnnp[k]) for n=1:lmax, k=1:nbK])` by default.
Using `StaticArrays` for better performances.
"""
const tabcosnf_Klnnp = SMatrix{lmax,nbK,Float64}([cos(n*tabf_Klnnp[k]) for n=1:lmax, k=1:nbK])

"""
    Table_K
    
`IntTable` structure used in the computation of the geometric factor Klnnp.
Use `Table_create!()` to initialize.
"""
const Table_K = Table_create!()

##################################################
# Computing the interaction coefficients Klnnp
##################################################
"""
    Klnnp(l,n,np,a,jt,ap,jp,Table_K_arg=Table_K)

Geometric factor K_nn'(a,j,a',j') involved in the computation of the scalar resonant diffusion coefficient DRR_jj.

# Remarks
    In practise this function does not exactly satisfy the symmetry Klnnp(a,j,ap,jp)=Klnnp(ap,jp,a,j).

# Arguments
- `l   ::Int64`  : mode of the multipole expansion.
- `n   ::Int64`  : first resonance number of the orbit.
- `np  ::Int64`  : second resonance number of the orbit.
- `a   ::Float64`: semi-major axis of the first orbit.
- `jt  ::Float64`: reduced angular momentum of the first orbit.
- `ap  ::Float64`: semi-major axis of the second orbit.
- `jp  ::Float64`: reduced angular momentum of the second orbit.
- `Table_K_arg=Table_K:: IntTable`: Structure which contains the tables needed for the computation of this integral.
"""
function Klnnp(l::Int64,n::Int64,np::Int64,a::Float64,j::Float64,ap::Float64,jp::Float64,Table_K_arg=Table_K)
    #####
    # Computing tabdeltaP and tabdeltaQ
    for j=1:nbK
        #####
        # Computing tabdeltaP
        x = 0.0
        for i=(Table_K_arg.tabw[j]+1):(Table_K_arg.tabw[j+1]) # ATTENTION, indices are shifted by one
            x += Table_K_arg.tabg[i]*(Table_K_arg.tabR[i]/Table_K_arg.tabRp[j])^(l)/(Table_K_arg.tabRp[j]) # Adding the contribution
        end
        Table_K_arg.tabdeltaP[j] = x # Updating tabdeltaP
        #####
        # Computing tabdeltaQ
        x = 0.0
        for i=(Table_K_arg.tabw[j+1]+1):(Table_K_arg.tabw[j+2]) # ATTENTION, indices are shifted by one
            x += Table_K_arg.tabg[i]*(Table_K_arg.tabRp[j]/Table_K_arg.tabR[i])^(l)/(Table_K_arg.tabR[i]) # Adding the contribution
        end
        Table_K_arg.tabdeltaQ[j] = x # Updating tabdeltaQ
    end
    #####
    # Computing the contributions from the sum P
    x = Table_K_arg.tabdeltaP[1] # Initial value of tabP[1]
    sumP = Table_K_arg.tabgp[1]*x # Initialising the contribution from P with the value of tabP[1]
    for j=1:(nbK-1) # Loop over the samples
        x = (Table_K_arg.tabRp[j]/Table_K_arg.tabRp[j+1])^(l+1)*x + Table_K_arg.tabdeltaP[j+1] # Computing the value of P[j+1]
        sumP += Table_K_arg.tabgp[j+1]*x # Updating the total contribution from P with P[j+1]
    end
    #####
    # Computing the contributions from the sum Q
    x = Table_K_arg.tabdeltaQ[nbK] # Initial value of tabQ[1]
    sumQ = Table_K_arg.tabgp[nbK]*x # Initialising the contribution from Q with the value of tabQ[nbK]
    for j=nbK:-1:2 # Loop over the samples. ATTENTION, this is a decreasing recurrence
        x = (Table_K_arg.tabRp[j-1]/Table_K_arg.tabRp[j])^(l)*x + Table_K_arg.tabdeltaQ[j-1] # Computing the value of Q[j-1]
        sumQ += Table_K_arg.tabgp[j-1]*x # Updating the total contribution from Q with Q[j-1]
    end
    #####
    return (sumP+sumQ)/(nbK^(2)) # Computing the total contribution from both P and Q
end



##################################################
# Computing the wrapped interaction coefficients
##################################################
"""
    AnnpSQ(l,n,np,a,jt,ap,jp,Table_K_arg=Table_K)

Coefficient |A_nnp|^2 at (a,jt,a',j') involved in the computation of the scalar resonant diffusion coefficient DRR_jj.

# Remark
    Do not forget that this function returns the SQUARE of the coefficient!

# Arguments
- `l   ::Int64`  : mode of the multipole expansion.
- `n   ::Int64`  : first resonance number of the orbit.
- `np  ::Int64`  : second resonance number of the orbit.
- `a   ::Float64`: semi-major axis of the first orbit.
- `jt  ::Float64`: reduced angular momentum of the first orbit.
- `ap  ::Float64`: semi-major axis of the second orbit.
- `jp  ::Float64`: reduced angular momentum of the second orbit.
- `Table_K_arg=Table_K:: IntTable`: Structure which contains the tables needed for the computation of this integral.
"""
function AnnpSQ(n::Int64,np::Int64,a::Float64,j::Float64,ap::Float64,jp::Float64,Table_K_arg=Table_K)
    # Getting rid of the cases where the coefficients is necessarily zero
    if (abs(n) > lmax || abs(np) > lmax || mod(n-np,2) != 0)
        return 0.0
    else # In that case the result is non-zero
        res = 0.0 # Initialising the result
        #####
        Table_init!(Table_K_arg,n,np,a,j,ap,jp) # Updating tabR, tabRp, tabg, tabgp, and tabw for the respective orders of the particles
        #####
        for l = max(abs(n),abs(np)):2:lmax # Loop over the l that have non-vanishing contributions !! ATTENTION, the step of the loop is 2
            res += (ylnSQ(l,n)*ylnSQ(l,np))/((2.0*l+1.0)^(3))*(Klnnp(l,n,np,a,j,ap,jp,Table_K_arg))^(2) # Adding the contribution from (l,n,np) !! ATTENTION, not to forget the squared
        end
        #####
        res *= 16.0*pi^(2) # Putting the prefactor
        return res # Output
    end
end
