TITLE Stochastic version of K-A channel from Klee Ficker and Heinemann 
: modified to account for Dax A Current --- M.Migliore Jun 1997
: modified to be used with cvode  M.Migliore 2001
: modified to be stochastic by K. Diba (changed original as little as possible)
  

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
 :   (pS) = (picosiemens)
 :   (um) = (micron)
}

PARAMETER {
    v       (mV)
    dt      (ms)
    area
    celsius     (degC)
    
    gamma =  10 (pS) : single channel conductance
    eta   =  8  (1/um2) : channel density 
    deterministic = 0   : if non-zero, use deterministic variables  
        
        vhalfn=11   (mV)
        vhalfl=-56   (mV)
        a0l=0.05      (/ms)
        a0n=0.05    (/ms)
        zetan=-1.5    (1)
        zetal=3    (1)
        gmn=0.55   (1)
        gml=1   (1)
    lmin=2  (mS)
    nmin=0.1  (mS)
    pw=-1    (1)
    tq=-40
    qq=5
    q10=5
    qtl=1
    ek                                    
    vmin = -120 (mV)    : range to construct tables for
    vmax = 100  (mV)
    DONT_VECTORIZE          : required declaration
}


NEURON {
    SUFFIX ska
    USEION k READ ek WRITE ik
    RANGE N, eta, gamma, gka, deterministic,reff
    GLOBAL P_an, P_bn, P_al, P_bl
    GLOBAL ninf, linf, ltau, ntau, lmin
    GLOBAL vmin, vmax, q10
    GLOBAL DONT_VECTORIZE   : prevent vectorization to agree with RNG.mod
}

STATE {
    n l  : deterministic state
    N0L0 N1L0 N0L1 N1L1 : states
    n0l0_n1l0 n0l0_n0l1 : transitions
    n1l0_n1l1 n1l0_n0l0
    n0l1_n1l1 n0l1_n0l0
    n1l1_n0l1 n1l1_n1l0
}

ASSIGNED {
    ik (mA/cm2)
        ninf
        linf      
        ltau  (ms)
        ntau   (ms)
        gka  (pS/um2)  
            an      (/ms)
    bn      (/ms)     al      (/ms)
    bl      (/ms)     reff    (pS/um2)

    N 
    scale_dens (pS/um2) 
    P_an     : probability of one channel making alpha n transition
    P_bn     : probability of one channel making beta n transition
    P_al     : probability of one channel making alpha l transition
    P_bl     : probability of one channel making beta l transition

}

INITIAL {
    rates(v)
    n=ninf
    l=linf
    scale_dens = gamma/area
    reff = eta*gamma
    N = floor(eta*area + 0.5)
    
    N1L1 = floor(n* l* N + 0.5)
    N1L0 = floor(n* (1-l)* N + 0.5)
    N0L1 = floor((1-n)* l* N + 0.5)
    N0L0 = N - N1L1 - N1L0 - N0L1  : put rest into non-conducting state
    
    n0l0_n1l0 = 0
    n0l0_n0l1 = 0
    n1l0_n1l1 = 0
    n1l0_n0l0 = 0
    n0l1_n1l1 = 0
    n0l1_n0l0 = 0
    n1l1_n0l1 = 0
    n1l1_n1l0 = 0
}


BREAKPOINT {
    SOLVE states 
    if (deterministic) { 
        if (deterministic-1){   
      gka = n * l * reff     
      } else {                                    
      gka = floor(n*l * N + 0.5) * scale_dens}
    } else{                                           
      gka = strap(N1L1) * scale_dens
    }
    ik = gka*(v-ek)*(1e-4)
}


: ----------------------------------------------------------------
: states - compute state variables
PROCEDURE states() {

VERBATIM
    extern double BnlDev_RNG();
ENDVERBATIM
        
    rates(v)

    : deterministic versions of state variables
    : integrated by relaxing toward the steady state value
    n = n + (1 - exp(-dt/ntau)) * (ninf-n)
    l = l + (1 - exp(-dt/ltau)) * (linf-l)
    
    P_an = strap(an*dt)
    P_bn = strap(bn*dt)
    
    : check that will represent probabilities when used
    ChkProb( P_an)
    ChkProb( P_bn)
    
    : n gate transitions
    
    n0l0_n1l0 = BnlDev_RNG(P_an, N0L0)    
    n0l1_n1l1 = BnlDev_RNG(P_an, N0L1)
    n1l1_n0l1 = BnlDev_RNG(P_bn, N1L1)
    n1l0_n0l0 = BnlDev_RNG(P_bn, N1L0)

    : new numbers in each state after the n gate transitions
    N0L0 = N0L0 - n0l0_n1l0 + n1l0_n0l0
    N1L0 = N1L0 - n1l0_n0l0 + n0l0_n1l0
    
    N0L1 = N0L1 - n0l1_n1l1 + n1l1_n0l1
    N1L1 = N1L1 - n1l1_n0l1 + n0l1_n1l1

    : probabilities of making l gate transitions
    P_al = strap(al*dt)
    P_bl  = strap(bl*dt)
    
    ChkProb(P_al)
    ChkProb(P_bl)
    
    : number making l gate transitions

    n0l0_n0l1 = BnlDev_RNG(P_al,N0L0-n0l0_n1l0)
    n1l0_n1l1 = BnlDev_RNG(P_al,N1L0-n1l0_n0l0)
    n0l1_n0l0 = BnlDev_RNG(P_bl,N0L1-n0l1_n1l1)
    n1l1_n1l0 = BnlDev_RNG(P_bl,N1L1-n1l1_n0l1)

    N0L0 = N0L0 - n0l0_n0l1  + n0l1_n0l0
    N1L0 = N1L0 - n1l0_n1l1  + n1l1_n1l0
    
    N0L1 = N0L1 - n0l1_n0l0  + n0l0_n0l1
    N1L1 = N1L1 - n1l1_n1l0  + n1l0_n1l1
 }

FUNCTION alpn(v(mV)) {
LOCAL zeta
  zeta=zetan+pw/(1+exp((v-tq)/qq))
  alpn = exp(1.e-3*zeta*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betn(v(mV)) {
LOCAL zeta
  zeta=zetan+pw/(1+exp((v-tq)/qq))
  betn = exp(1.e-3*zeta*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION alpl(v(mV)) {
  alpl = exp(1.e-3*zetal*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betl(v(mV)) {
  betl = exp(1.e-3*zetal*gml*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius))) 
}    
PROCEDURE rates(v(mV)) { :callable from hoc
        LOCAL a,qt    
        TABLE ntau, ltau, ninf, linf,al,bl,an,bn
        DEPEND q10, celsius, a0n, nmin, lmin,qtl
        FROM vmin TO vmax WITH 199
         
        qt=q10^((celsius-24)/10)
        a = alpn(v)
        ninf = 1/(1 + a)
        ntau = betn(v)/(qt*a0n*(1+a))
    if (ntau<nmin) {ntau=nmin}
        a = alpl(v)
        linf = 1/(1+ a)
    ltau = 0.26*(v+50)/qtl
    if (ltau<lmin/qtl) {ltau=lmin/qtl}   
    al = linf/ltau
    bl = 1/ltau - al
    an = ninf/ntau
    bn = 1/ntau - an
}

: ----------------------------------------------------------------
: sign trap - trap for negative values and replace with zero
FUNCTION strap(x) {
    if (x < 0) {
        strap = 0
VERBATIM
        fprintf (stderr,"skaprox.mod:strap: negative state");
ENDVERBATIM
    } else {
        strap = x
    }
}

: ----------------------------------------------------------------
: ChkProb - Check that number represents a probability
PROCEDURE ChkProb(p) {
  if (p < 0.0 || p > 1.0) {
VERBATIM
    fprintf(stderr, "skaprox.mod:ChkProb: argument not a probability.\n");
ENDVERBATIM
  }
}


