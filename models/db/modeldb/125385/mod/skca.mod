TITLE skca.mod  
 
COMMENT
----------------------------------------------------------------
Stochastic version of kca.mod

Calcium-dependent potassium channel
Based on
Pennefather (1990) -- sympathetic ganglion cells
taken from
Reuveni et al (1993) -- neocortical cells

Author: Zach Mainen, Salk Institute, 1995, zach@salk.edu
                   
19 May 2002 Kamran Diba.  Changed gamma and deterministic from GLOBAL to RANGE.
----------------------------------------------------------------
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX skca
    USEION k READ ek WRITE ik
    USEION ca READ cai        
    RANGE gk, gamma, eta, N, deterministic,reff
    GLOBAL P_a,P_b, ninf, ntau,a,b     
    GLOBAL Ra, Rb, caix
    GLOBAL vmin, vmax, q10, temp, tadj
    GLOBAL DONT_VECTORIZE   : prevent vectorization to agree with RNG.mod
}
   
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (um) = (micron)
} 

PARAMETER {
    v           (mV)
    dt      (ms)
    area
    
    gamma = 180       (pS)    : 0.03 mho/cm2
    eta = 0.0556         (1/um2)     
    cai         (mM)
    caix = 1    
                                    
    Ra   = 0.01 (/ms)       : max act rate  
    Rb   = 0.02 (/ms)       : max deact rate 
        celsius     (degC)
    temp = 23   (degC)      : original temp     
    q10  = 2.3          : temperature sensitivity
    deterministic = 0   : if non-zero, will use deterministic version
    vmin = -120 (mV)    : range to construct tables for
    vmax = 100  (mV)
    DONT_VECTORIZE      : required declaration
} 

ASSIGNED {
    a       (/ms)
    b       (/ms)
    ik      (mA/cm2)
    gk      (pS/um2)
    ek      (mV)
    ninf        : steady-state value
    ntau (ms)   : time constant for relaxation    
    tadj              

    N 
    reff    (pS/um2)
    scale_dens (pS/um2) 
    P_a     : probability of one channel making alpha transition
    P_b     : probability of one channel making beta transition
}
 

STATE {
    n         : state variable of deterministic description
    N0 N1       : N states populations
    n0_n1 n1_n0 : number of channels moving from one state to the other 
}
: ----------------------------------------------------------------
: initialization
INITIAL { 
    rates(cai)
    n = ninf
    scale_dens = gamma/area
    reff = eta*gamma
    N = floor(eta*area + 0.5)
    
    N1 = floor(n * N + 0.5)
    N0 = N-N1       : any round off into non-conducting state
    
    n0_n1 = 0
    n1_n0 = 0
}

: ----------------------------------------------------------------
: Breakpoint for each integration step
BREAKPOINT {
  SOLVE states
   if (deterministic) { 
        if (deterministic-1){
    gk =  n *reff * tadj
    } else { 
    gk = floor(n* N + 0.5) * scale_dens *tadj} 
    } else{                                         
    gk =  strap(N1) * scale_dens * tadj
    }
    ik = (1e-4) * gk * (v - ek)    
}


: ----------------------------------------------------------------
: states - updates number of channels in each state
PROCEDURE states() {

VERBATIM
    extern double BnlDev_RNG();
ENDVERBATIM

    rates(cai)
    
    : deterministic versions of state variables
    : integrated by relaxing toward steady-state
    n = n + (1 - exp(-dt/ntau)) * (ninf-n)

    P_a = strap(a*dt)
    P_b = strap(b*dt)

    : check that will represent probabilities when used
    ChkProb( P_a)
    ChkProb( P_b)
    
    : transitions
:    if (deterministic) {
:    n0_n1 = P_a*N0
:    n1_n0 = P_b*N1
:    }
:    else{    
    n0_n1 = BnlDev_RNG(P_a, N0)
    n1_n0 = BnlDev_RNG(P_b, N1)
:    }
    : move the channels
    N0    = strap(N0 - n0_n1 + n1_n0)
    N1    = N - N0
}

PROCEDURE rates(cai(mM)) {
        
        tadj = q10^((celsius - temp)/10)
                                  
        a = Ra * cai^caix
    a = a * tadj                                     
        b = Rb
    b = b * tadj                                      
        ntau = 1/(a+b)
        ninf = a/(a+b)
}                                   

: ----------------------------------------------------------------
: sign trap - trap for negative values and replace with zero
FUNCTION strap(x) {
    if (x < 0) {
        strap = 0
VERBATIM
        fprintf (stderr,"skca.mod:strap: negative state");
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
    fprintf(stderr, "skca.mod:ChkProb: argument not a probability.\n");
ENDVERBATIM
  }
}
