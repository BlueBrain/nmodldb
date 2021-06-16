TITLE skm95.mod  
 
COMMENT
----------------------------------------------------------------
Stochastic version of km.mod

Potassium channel, Hodgkin-Huxley style kinetics
Based on I-M (muscarinic K channel)
Slow, noninactivating

Author: Zach Mainen, Salk Institute, 1995, zach@salk.edu
19 May 2002 Kamran Diba.  Changed gamma and deterministic from GLOBAL to RANGE.
----------------------------------------------------------------
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (um) = (micron)
} 

NEURON {
    SUFFIX skm
    USEION k READ ek WRITE ik
    RANGE gk, gamma, deterministic,reff
    RANGE N,eta
    GLOBAL ninf, ntau,a,b, P_a,P_b
    GLOBAL Rb, Ra
    GLOBAL vmin, vmax, q10, temp, tadj
    GLOBAL DONT_VECTORIZE   : prevent vectorization to agree with RNG.mod
}
   
PARAMETER {
    v           (mV)
    dt      (ms)
    area
    
    gamma = 40       (pS)    : 0.03 mho/cm2
    eta = 0.25       (1/um2)                                 
    tha  = -30  (mV)        : v 1/2 for inf
    qa   = 9    (mV)        : inf slope       
    
    Ra   = 0.001    (/ms)       : max act rate  (slow)
    Rb   = 0.001    (/ms)       : max deact rate  (slow)
    
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
    trates(v)
    n = ninf
    scale_dens = gamma/area
    reff = eta*gamma
    N = floor(eta*area + 0.5)
        
    N1 = floor(n * N+ 0.5)
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
    gk = n*reff*tadj
    } else {
    gk = floor(n*N + 0.5) * scale_dens*tadj}
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

    trates(v)
    
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

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

        TABLE ninf, ntau, a, b, tadj
    DEPEND  celsius, temp, Ra, Rb, tha, qa,q10
    FROM vmin TO vmax WITH 199

        tadj = q10^((celsius - temp)/10) 
        a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa))
    a = a * tadj                    
        b = -Rb * (v - tha) / (1 - exp((v - tha)/qa))
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
        fprintf (stderr,"skm.mod:strap: negative state");
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
    fprintf(stderr, "skm.mod:ChkProb: argument not a probability.\n");
ENDVERBATIM
  }
}
