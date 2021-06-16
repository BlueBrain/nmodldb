TITLE nax
: Na current for axon. No slow inact.
: M.Migliore Jul. 1997
: added sh to account for higher threshold M.Migliore, Apr.2002
: MM: in parallel to updated na.mod from Mainen & Sejnowski,1996
: trap0 was replaced with efun, which is a better approximation
: in the vicinity of a singularity
: MM, mtaufac added to explore the effect of faster activation
: MM, htaufac added to model faster inactivation of axonal Na-channels

NEURON {
    SUFFIX nax
    USEION na READ ena WRITE ina
    RANGE  gbar, sh, shx, mtaufac, htaufac, qa, inax, thegna, m, h, m3h, qinf, thinf
    GLOBAL minf, hinf, mtau, htau
}

PARAMETER {
    sh   = -7   (mV)        : leftward shift of minf and mtau
    shx  = -7    (mV)        : shift of hinf and htau
    gbar = 0.010    (mho/cm2)   
    mtaufac = 1     (1)     : factor that is multiplied with the expression for mtau    
    htaufac = 1     (1)     : factor that is multiplied with the expression for htau    
                            
    tha  =  -30 (mV)        : v 1/2 for act 
    qa   = 7.2  (mV)        : act slope (4.5)       
    Ra   = 0.4  (/ms)       : open (v)      
    Rb   = 0.124    (/ms)       : close (v)     

    thi1  = -45 (mV)        : v 1/2 for inact   
    thi2  = -45     (mV)        : v 1/2 for inact   
    qd   = 1.5  (mV)            : inact tau slope
    qg   = 1.5      (mV)
    mmin=0.02   
    hmin=0.5            
    q10=2
    Rg   = 0.01     (/ms)       : inact recov (v)   
    Rd   = .03  (/ms)       : inact (v) 

    thinf  = -50    (mV)        : inact inf V_1/2 - originally -50 (changed to -62)
    qinf  =  4  (mV)        : inact inf slope - originally 4 (changed to 6.9)

    ena     (mV)            : must be explicitly def. in hoc
    celsius 
    v       (mV)
}


UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (um) = (micron)
} 

ASSIGNED {
    ina         (mA/cm2)
    inax        (mA/cm2)
    thegna      (mho/cm2)
    minf        hinf
    m3h
    mtau (ms)   htau (ms)   
}
 

STATE { m h}

BREAKPOINT {
        SOLVE states METHOD cnexp
        thegna = gbar*m*m*m*h
        m3h = m*m*m*h
    ina = thegna * (v - ena)
    inax = thegna * (v - ena)
    ina = inax
} 

INITIAL {
    trates(v,sh, shx, mtaufac, htaufac, qa, qinf)
    m=minf  
    h=hinf
}

DERIVATIVE states {   
        trates(v,sh,shx, mtaufac, htaufac, qa, qinf)      
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
}

: efun() is a better approx than trap0 in vicinity of singularity--

PROCEDURE trates(vm,sh2, sh3, taufac, taufac2, qa, qinf) {  
        LOCAL  a, b, qt
        qt=q10^((celsius-24)/10)
:   a = trap0(vm,tha+sh2,Ra,qa)
    a = Ra * qa * efun((tha+sh2 - vm)/qa)
:   b = trap0(-vm,-tha-sh2,Rb,qa)
    b = Rb * qa * efun((vm - tha-sh2)/qa)

:    mtau = taufac*0.001    
:   mtau = 1/(a+b)/qt
:        if (mtau<mmin) {mtau=mmin}
    mtau = taufac/(a+b)/qt
        if (mtau<mmin) {mtau=taufac*mmin}
    minf = a/(a+b)
    
:   a = trap0(vm,thi1+sh2,Rd,qd)
    a = Rd * qd * efun((thi1+sh3 - vm)/qd)
:   b = trap0(-vm,-thi2-sh2,Rg,qg)
    b = Rg * qg * efun((vm - thi2-sh3)/qg)
    
:   htau =  1/(a+b)/qt
:        if (htau<hmin) {htau=hmin}
    htau =  taufac2/(a+b)/qt
        if (htau<hmin) {htau=taufac2*hmin}
    hinf = 1/(1+exp((vm-thinf-sh3)/qinf))
}
        
COMMENT
FUNCTION trap0(v,th,a,q) {
    if (fabs(v-th) > 1e-6) {
            trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
    } else {
            trap0 = a * q
    }
}   
ENDCOMMENT

FUNCTION efun(z) {
    if (fabs(z) < 1e-6) {
        efun = 1 - z/2
    }else{
        efun = z/(exp(z) - 1)
    }
}        
