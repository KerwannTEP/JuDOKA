##################################################
# Computation of the NR diffusion coefficients
##################################################

##################################################
# Preliminary useful functions
##################################################

"""
    logCoulomb(iBath)
   
Coulomb logarithm of the family indexed by iBath in the source file.
Defined as the logarithm of the ratio of the mass of the SMBH over the mass of an element of the family iBath.

# Remark
Using the convention where log is the Neperian logarithm (in base e).

"""
function logCoulomb(iBath::Int64)
    m = tabmBath[iBath]
    Q = mBH/m
    return log(Q)
end

"""
    prefac_J(E)
   
Numerical prefactor J(E) used to convert the PDF f(E,R) to the number density N(E,R), where R = j^2.
Depends on the binding energy E, which is the opposite of the gravitational energy.

Explicit formula is N(E,R) = J(E)*f(E,R).
"""
function prefac_J(E::Float64)
    return (2.0)^(1/2)*pi^3*(G*mBH)^3*E^(-5/2)
end

"""
    E_bind(a)
    
Binding energy of the system, which is the opposite of the gravitational energy.
Defined as `E(a) = G*mBH / (2a)`, where `a` is the semi-major axis of the orbit, G is the gravitational constant and mBH the mass of the SMBH.
"""
function E_bind(a::Float64)
    return (G*mBH)/(2.0*a)
end

"""
    j_to_R(j)
    
Converts reduced angular momentum j = J/Jc(a) into R = j^2, where J is the angular momentum and Jc(a) the angular momentum of a circular orbit.
"""
function j_to_R(j::Float64) 
    return j^2
end

"""
    DF_E(E,iBath)
    
R-averaged distribution function in binding enery E, defined in E space, for the family indexed by iBath in the source file.
R is the squared reduced angular momentum j^2.
"""
function DF_E(E::Float64,iBath::Int64)
    E0    = E_bind(taba0Bath[iBath])
    NBar  = tabNbarBath[iBath]
    gamma = tabgammaBath[iBath]
    return NBar/prefac_J(E) * (G*mBH)/(2.0*E^2) * (E0/E)^(2.0-gamma)
end

"""
    gamma_prefactor(iBath)
    
Gamma prefactor which appears in Fi integrals, for the family indexed by iBath in the source file.
Used in the computation of the non-resonant diffusion coefficients.
"""
function gamma_prefactor(iBath::Int64)
    m = tabmBath[iBath]
    return 4.0*pi*(G*m)^2 * logCoulomb(iBath)
end

"""
    xpm(R)

xplus and xminus quantities which appear in the bounds of the integrals involved in the computation of the non-resonant diffusion coefficients. 

# Remark
x = E(r)/gravitational potential(r)
"""
function xpm(R::Float64)
    return (1.0+sqrt(1.0-R))/2.0,(1.0-sqrt(1.0-R))/2.0
end

"""
    inverse_reduced_E_bind(Ep,E)

Reduced binding energy used in the computation of integrals Fi.
Defined as Ep over E.
"""
function inverse_reduced_E_bind(Ep::Float64,E::Float64)
    return Ep/E
end

##################################################
# Bounds in A_k, B_k and C_k (xp>1/s) integrals
##################################################

"""
    _beta(R,s)

Bound in Ck integral when x_plus > 1/s, where s is the reduced binding energy and R is the squared reduced angular momentum.
"""
function _beta(R::Float64,s::Float64)
    xp,xm = xpm(R)
    res = 1.0/(xm*s)
    res = sqrt(res)
    return acosh(res)
end

"""
    _eta(R)

Bound in Ak integral when xp <= 1/s, where s is the reduced binding energy and R is the squared reduced angular momentum.
"""
function _eta(R::Float64)
    xp,xm = xpm(R)
    res = 1.0/(2.0*xm)
    res = sqrt(res)
    return acosh(res)
end

"""
    _delta(R)

Bound in Bk integral when xp <= 1/s, where s is the reduced binding energy and R is the squared reduced angular momentum.
"""
function _delta(R::Float64)
    xp,xm = xpm(R)
    res = 1.0/(2.0*xp)
    res = sqrt(res)
    return acos(res)
end

##################################################
# Computation of the NR diffusion coefficients.
##################################################

"""
    _fB(theta,s,R)

Computes all seven integrands within Bk integrals for 1<=k<=7.
Returns a 7-uple, ordered from k=1 to k=7.
See `Conrad et Kulsrud (1978)`.

# Remarks
- All the integrands are computed at the same time.
- Avoids multiple evalution of cos(theta), xplus and xminus.
- Useful when xplus <= 1/s.

# Arguments
- `theta::Float64`: Integration variable.
- `s    ::Float64`: Reduced binding energy.
- `R    ::Float64`: Squared reduced angular momentum j^2.
"""
function _fB(theta::Float64,s::Float64,R::Float64)
    xp,xm = xpm(R)
    ym    = xm/xp
    csTh  = cos(theta)

    fB1 = csTh^3/sqrt(csTh^2-ym) * sqrt(1.0-s*xp*csTh^2)^(1)/sqrt(1.0-xp*csTh^2)^(1)
    fB2 = csTh^5/sqrt(csTh^2-ym) * sqrt(1.0-s*xp*csTh^2)^(1)/sqrt(1.0-xp*csTh^2)^(3)
    fB3 = csTh^7/sqrt(csTh^2-ym) * sqrt(1.0-s*xp*csTh^2)^(1)/sqrt(1.0-xp*csTh^2)^(1)
    fB4 = csTh/sqrt(csTh^2-ym)   * sqrt(1.0-s*xp*csTh^2)^(3)/sqrt(1.0-xp*csTh^2)^(1)
    fB5 = csTh^3/sqrt(csTh^2-ym) * sqrt(1.0-s*xp*csTh^2)^(3)/sqrt(1.0-xp*csTh^2)^(3)
    fB6 = csTh^5/sqrt(csTh^2-ym) * sqrt(1.0-s*xp*csTh^2)^(3)/sqrt(1.0-xp*csTh^2)^(5)
    fB7 = csTh^7/sqrt(csTh^2-ym) * sqrt(1.0-s*xp*csTh^2)^(3)/sqrt(1.0-xp*csTh^2)^(3)
    return fB1,fB2,fB3,fB4,fB5,fB6,fB7
end

"""
    _fAC(theta,s,R)

Computes all 14 integrands within Ak and Ck integrals for 1<=k<=7.
Returns a 7-uple, ordered from k=1 to k=7.
See `Conrad et Kulsrud (1978)`.

# Remarks
- All the integrands are computed at the same time.
- Avoids multiple evalution of cos(theta), xplus and xminus.
- Useful when xplus <= 1/s for Ak.
- Useful when xplus > 1/s for Ck.

# Arguments
- `theta::Float64`: Integration variable.
- `s    ::Float64`: Reduced binding energy.
- `R    ::Float64`: Squared reduced angular momentum j^2.
"""
function _fAC(theta::Float64,s::Float64,R::Float64)
    xp,xm = xpm(R)
    ym    = xp/xm
    cshTh = cosh(theta)

    fAC1 = cshTh^3/sqrt(ym-cshTh^2) * sqrt(1.0-s*xm*cshTh^2)^(1)/sqrt(1.0-xm*cshTh^2)^(1)
    fAC2 = cshTh^5/sqrt(ym-cshTh^2) * sqrt(1.0-s*xm*cshTh^2)^(1)/sqrt(1.0-xm*cshTh^2)^(3)
    fAC3 = cshTh^7/sqrt(ym-cshTh^2) * sqrt(1.0-s*xm*cshTh^2)^(1)/sqrt(1.0-xm*cshTh^2)^(1)
    fAC4 = cshTh  /sqrt(ym-cshTh^2) * sqrt(1.0-s*xm*cshTh^2)^(3)/sqrt(1.0-xm*cshTh^2)^(1)
    fAC5 = cshTh^3/sqrt(ym-cshTh^2) * sqrt(1.0-s*xm*cshTh^2)^(3)/sqrt(1.0-xm*cshTh^2)^(3)
    fAC6 = cshTh^5/sqrt(ym-cshTh^2) * sqrt(1.0-s*xm*cshTh^2)^(3)/sqrt(1.0-xm*cshTh^2)^(5)
    fAC7 = cshTh^7/sqrt(ym-cshTh^2) * sqrt(1.0-s*xm*cshTh^2)^(3)/sqrt(1.0-xm*cshTh^2)^(3)
    return fAC1,fAC2,fAC3,fAC4,fAC5,fAC6,fAC7
end

"""
    _C_NR(s,R,nbC=500) 

Computes all seven Ck integrals for 1<=k<=7 used for the computation of the non-resonance diffusion coefficients.
Midpoint sampling rule is used for integration.
Returns a 7-uple, ordered from k=1 to k=7.
See `Conrad et Kulsrud (1978)`.

# Remarks
- All the integrands are computed at the same time.
- Avoids multiple evalution of cos(theta), xplus and xminus, etc.

# Arguments
- `s  ::Float64`: Reduced binding energy.
- `R  ::Float64`: Squared reduced angular momentum j^2.
- `nbC::Float64`: Sampling number of the computed integrals. Set as 500 by default.
"""
function _C_NR(s::Float64,R::Float64,nbC::Int64=500) 
    xp,xm = xpm(R)
    C1,C2,C3,C4,C5,C6,C7 = 0.0,0.0,0.0,0.0,0.0,0.0,0.0

    if xp>1/s
        supC = _beta(R,s)
        stepC = supC/nbC
        prefactorC = 4.0*stepC/pi
    else
        supA = _eta(R)
        supB = _delta(R)
        stepA = supA/nbC
        stepB = supB/nbC
        prefactorA = 4.0*stepA/pi
        prefactorB = 4.0*stepB/pi
        A1,A2,A3,A4,A5,A6,A7 = 0.0,0.0,0.0,0.0,0.0,0.0,0.0
        B1,B2,B3,B4,B5,B6,B7 = 0.0,0.0,0.0,0.0,0.0,0.0,0.0
    end
    for jTheta=1:nbC
        if xp>1/s
            thetaC = (jTheta-0.5)*stepC
            fAC1,fAC2,fAC3,fAC4,fAC5,fAC6,fAC7 = _fAC(thetaC,s,R)
            C1 += fAC1
            C2 += fAC2
            C3 += fAC3
            C4 += fAC4
            C5 += fAC5
            C6 += fAC6
            C7 += fAC7
        else
            thetaA = (jTheta-0.5)*stepA
            thetaB = (jTheta-0.5)*stepB
            fAC1,fAC2,fAC3,fAC4,fAC5,fAC6,fAC7 = _fAC(thetaA,s,R)
            fB1,fB2,fB3,fB4,fB5,fB6,fB7 = _fB(thetaB,s,R)
            A1 += fAC1
            A2 += fAC2
            A3 += fAC3
            A4 += fAC4
            A5 += fAC5
            A6 += fAC6
            A7 += fAC7
            B1 += fB1
            B2 += fB2
            B3 += fB3
            B4 += fB4
            B5 += fB5
            B6 += fB6
            B7 += fB7
        end
    end
    if xp>1/s
        C1 *= prefactorC*xm
        C2 *= prefactorC*xm^2
        C3 *= prefactorC*xm^3
        C4 *= prefactorC
        C5 *= prefactorC*xm
        C6 *= prefactorC*xm^2
        C7 *= prefactorC*xm^3
    else
        A1 *= prefactorA*xm
        A2 *= prefactorA*xm^2
        A3 *= prefactorA*xm^3
        A4 *= prefactorA
        A5 *= prefactorA*xm
        A6 *= prefactorA*xm^2
        A7 *= prefactorA*xm^3
        B1 *= prefactorB*xp
        B2 *= prefactorB*xp^2
        B3 *= prefactorB*xp^3
        B4 *= prefactorB
        B5 *= prefactorB*xp
        B6 *= prefactorB*xp^2
        B7 *= prefactorB*xp^3
        C1 = A1 + B1
        C2 = A2 + B2
        C3 = A3 + B3
        C4 = A4 + B4
        C5 = A5 + B5
        C6 = A6 + B6
        C7 = A7 + B7
    end
    return C1,C2,C3,C4,C5,C6,C7
end

"""
    _F_NR_0(E,iBath)
    
Analytic evaluation of the F0 integral, at binding energy E and for the family indexed by iBath in the source file.
See `Conrad et Kulsrud (1978)`.
"""
function _F_NR_0(E::Float64,iBath::Int64)
    a0 = taba0Bath[iBath]
    E0 = E_bind(a0)
    NBar  = tabNbarBath[iBath]
    gamma = tabgammaBath[iBath]
    gamma_pref = gamma_prefactor(iBath)
    return 4.0*pi*gamma_pref* a0 * 1.0/prefac_J(E)*NBar/(gamma-0.5) * (E0/E)^(3.0-gamma)
end

"""
    _F_NR(E,R,iBath,nbC=500) 

Computes all seven Fk integrals for 1<=k<=7 used for the computation of the non-resonance diffusion coefficients.
Midpoint sampling rule is used for integration.
Returns a 7-uple, ordered from k=1 to k=7.
See `Conrad et Kulsrud (1978)`.

# Remarks
- All the integrands are computed at the same time.
- Avoids multiple evalution of cos(theta), xplus and xminus, etc.
- Logarithmic sampling in energy for better precision at low values of R.

# Arguments
- `E    ::Float64`: Binding energy.
- `R    ::Float64`: Squared reduced angular momentum j^2.
- `iBath::Float64`: Index in the source file of the considered family.
- `nbC  ::Float64`: Sampling number of the computed integrals. Set as 500 by default.
"""
function _F_NR(E::Float64,R::Float64,iBath::Int64,nbF::Int64=500) 
    xp,xm = xpm(R)
    inf = E
    sup = E/xm
    log_step = log(sup / inf) / nbF
    # step = E/nbF * (1.0/xm - 1.0)
    F1,F2,F3,F4,F5,F6,F7 = 0.0,0.0,0.0,0.0,0.0,0.0,0.0
    for kEnergy=1:nbF
        Ep = inf * exp((kEnergy-0.5)*log_step)
        # Ep = inf+(kEnergy-0.5)*step
        s  = inverse_reduced_E_bind(Ep,E)
        C1,C2,C3,C4,C5,C6,C7 = _C_NR(s,R)
        DF = DF_E(Ep,iBath) * Ep
        F1 += DF*C1
        F2 += DF*C2
        F3 += DF*C3
        F4 += DF*C4
        F5 += DF*C5
        F6 += DF*C6
        F7 += DF*C7
    end
    prefactor = 4.0*pi*gamma_prefactor(iBath)*log_step
    # prefactor = 4.0*pi*gamma_prefactor(iBath)*step
    F1 *= prefactor
    F2 *= prefactor
    F3 *= prefactor
    F4 *= prefactor
    F5 *= prefactor
    F6 *= prefactor
    F7 *= prefactor
    return F1,F2,F3,F4,F5,F6,F7
end

##################################################
# Diffusion coefficients
##################################################

"""
    DNR_EJ(E,j,iBath,m_test=0.0)  

Non-resonant diffusion coefficients in (E,R) coordinates and Bar-Or&Alexander (2016) convention.
Returns the 5-uple `(DE, DEE/E, DJ/Jc, DJJ/Jc^2, DEJ/Jc)` for the family indexed by iBath.

# Arguments
- `E     ::Float64`: Binding energy.
- `j     ::Float64`: Reduced angular momentum `j=J/Jc`.
- `iBath ::Float64`: Index in the source file of the considered family.
- `m_test::Float64`: Mass of a test particle.
"""
function DNR_EJ(E::Float64,j::Float64,iBath::Int64,m_test::Float64=0.0)
    R = j^2
    F1,F2,F3,F4,F5,F6,F7 = _F_NR(E,R,iBath)
    F0 = _F_NR_0(E,iBath)
    m = tabmBath[iBath]
    D_E = -F0 + m_test/m * F1
    D_EE_over_E = 4/3*(F0 + F4)
    D_J_over_Jc = 1.0/(j*E) * ( (5.0 - 3.0*j^2)/(12.0)*F0 - j^2*(m_test + m)/(2.0*m)*F2 + F3 - 1/3*F7 )
    D_JJ_over_JcSQ = 1.0/E* ( (5.0 - 3.0*j^2)/6.0*F0 + 1/2*j^2*F6 - 1/2*j^2*F2 + 2.0*F3 - 2/3*F7 )
    D_EJ_over_Jc = -2/3*j*(F0 + F5)
    return D_E,D_EE_over_E,D_J_over_Jc,D_JJ_over_JcSQ,D_EJ_over_Jc
end

"""
    DNR_J_red(a,j)  

Total non-resonant reduced diffusion coefficient DJ/Jc in (E,R) coordinates and Bar-Or&Alexander (2016) convention.

# Arguments
- `a::Float64`: Semi-major axis.
- `j::Float64`: Reduced angular momentum `j=J/Jc`.
"""
function DNR_J_red(a::Float64,j::Float64)
    res = 0.0
    E = E_bind(a)
    for iBath=1:nbBath
        res += DNR_EJ(E,j,iBath)[3]
    end
    return res
end

"""
    DNR_JJ_red(a,j)  

Total non-resonant reduced diffusion coefficient DJJ/Jc^2 in (E,R) coordinates and Bar-Or&Alexander (2016) convention.

# Arguments
- `a::Float64`: Semi-major axis.
- `j::Float64`: Reduced angular momentum `j=J/Jc`.
"""
function DNR_JJ_red(a::Float64,j::Float64)
    res = 0.0
    E = E_bind(a)
    for iBath=1:nbBath
        res += DNR_EJ(E,j,iBath)[4]
    end
    return res
end

"""
    DNR_red(a,j,iBath)  

Non-resonant reduced diffusion coefficients in (a,j) coordinates and Bar-Or&Alexander (2016) convention.
Returns the 5-uple `(Da, Daa, Dj, Djj, Daj)` for the family indexed by iBath.

# Arguments
- `a     ::Float64`: Semi-major axis.
- `j     ::Float64`: Reduced angular momentum `j=J/Jc`.
- `iBath ::Float64`: Index in the source file of the considered family.
"""
function DNR_red(a::Float64,j::Float64,iBath::Int64) 
    E = E_bind(a)
    cst = 2.0/(G*mBH)
    D_E,D_EE_over_E,D_J_over_Jc,D_JJ_over_JcSQ,D_EJ_over_Jc = DNR_EJ(E,j,iBath)
    D_a = cst*a^2*(-D_E + D_EE_over_E)
    D_aa = cst*a^3*D_EE_over_E
    D_j = D_J_over_Jc + j/(2.0*E)*D_E + 1.0/(2.0*E)*D_EJ_over_Jc - j/(8.0*E)*D_EE_over_E
    D_jj = D_JJ_over_JcSQ + j/E*D_EJ_over_Jc + j^2/(4.0*E)*D_EE_over_E
    D_aj = -cst*a^2*(D_EJ_over_Jc + j/2*D_EE_over_E)
    return D_a,D_aa,D_j,D_jj,D_aj
end

"""
    DNR_aj(a,j)  

Total non-resonant reduced diffusion coefficients in (a,j) coordinates and Bar-Or&Alexander (2016) convention.
Returns the 5-uple `(Da, Daa, Dj, Djj, Daj)`.

# Arguments
- `a::Float64`: Semi-major axis.
- `j::Float64`: Reduced angular momentum `j=J/Jc`.
"""
function DNR_aj(a::Float64,j::Float64) #Wrapped function to give total D_NR_jj
    D_a,D_aa,D_j,D_jj,D_aj = 0.0,0.0,0.0,0.0,0.0
    for iBath=1:nbBath
        Di_a,Di_aa,Di_j,Di_jj,Di_aj = DNR_red(a,j,iBath)
        D_a  += Di_a
        D_aa += Di_aa
        D_j  += Di_j
        D_jj += Di_jj
        D_aj += Di_aj
    end
    return D_a,D_aa,D_j,D_jj,D_aj
end

"""
    DNR_jj(a,j)  

Total non-resonant reduced diffusion coefficient Djj in (a,j) coordinates and Bar-Or&Alexander (2016) convention.

# Arguments
- `a::Float64`: Semi-major axis.
- `j::Float64`: Reduced angular momentum `j=J/Jc`.
"""
function DNR_jj(a::Float64,j::Float64)
    return DNR_aj(a,j)[4]
end

"""
    DNR_j(a,j)  

Total non-resonant reduced diffusion coefficient Dj in (a,j) coordinates and Bar-Or&Alexander (2016) convention.

# Arguments
- `a::Float64`: Semi-major axis.
- `j::Float64`: Reduced angular momentum `j=J/Jc`.
"""
function DNR_j(a::Float64,j::Float64) 
    return DNR_aj(a,j)[3]
end