TITLE skm95.mod  
 
COMMENT
----------------------------------------------------------------
Stochastic version of the K channel mechanism kd3h5.mod by
Z. Mainen in Mainen & Sejnowski 95.

This represents a potassium channel, with Hodgkin-Huxley like kinetics,
based on the gates model, assuming stochastic opening and closing.

Kinetic rates based roughly on Sah et al. and Hamill et al. (1991)
The main kinetic difference from the standard H-H model (shh.mod) is 
that the K+ kinetic is different, not n^4, but just n, 
and the activation curves are different.

The rate functions are adapted directly from the Kd3h5.mod file
by Zach Mainen.

The stochastic model is as following:

Potassium

       = alpha_n =>      
   [N0]             [N1]
      <= beta_n =      


The model keeps track on the number of channels in each state, and 
uses a binomial distribution to update these number. This mechanism 
assumes that the RNG mechanism is also linked to provide binomial 
random deviates.

Jan 1999, Mickey London, Hebrew University, mikilon@lobster.ls.huji.ac.il
        Peter N. Steinmetz, Caltech, peter@klab.caltech.edu
14 Sep 99 PNS. Added deterministic flag.
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

PARAMETER {
    v           (mV)
    dt      (ms)
    area
    
    gamma  =  15      (pS)    
    eta      = 0.667      (1/um2)
    
    tha  = 25   (mV)        : v 1/2 for inf
    qa   = 9    (mV)        : inf slope     
    Ra   = 0.02 (/ms)       : max act rate
    Rb   = 0.002    (/ms)       : max deact rate
    
    celsius (degC)
    temp = 23 (degC)   : original temperature for kinetic set
    q10 = 2.3               : temperature sensitivity
    
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
    reff    (pS/um2)

    N 
    scale_dens (pS/um2) 
    P_a     : probability of one channel making alpha transition
    P_b     : probability of one channel making beta transition
}
 

STATE {
    n         : state variable of deterministic description
    N0 N1       : N states populations
    n0_n1 n1_n0 : number of channels moving from one state to the other 
}

NEURON {
    SUFFIX sk
    USEION k READ ek WRITE ik
    RANGE N,eta, gk, gamma, deterministic, reff     
    GLOBAL ninf, ntau,a,b,P_a,P_b
    GLOBAL Ra, Rb
    GLOBAL vmin, vmax, q10, temp, tadj
    GLOBAL DONT_VECTORIZE   : prevent vectorization to agree with RNG.mod
}

: ----------------------------------------------------------------
: initialization
INITIAL { 
    trates(v)
    n = ninf
    scale_dens = gamma/area
    N = floor(eta*area + 0.5)
    reff = eta*gamma
    
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

: ----------------------------------------------------------------
: trates - compute rates, using table if possible
PROCEDURE trates(v) {     
    TABLE ntau, ninf, a, b, tadj
    DEPEND dt, Ra, Rb, tha, qa, q10, temp, celsius
    FROM vmin TO vmax WITH 199
    
    tadj = q10 ^ ((celsius - temp)/10)
    a = SigmoidRate(v, tha, Ra, qa)
    a = a * tadj
    b = SigmoidRate(-v, -tha, Rb, qa)
    b = b * tadj
    ntau = 1/(a+b)
    ninf = a*ntau
}


: ----------------------------------------------------------------
: SigmoidRate - Compute a sigmoid rate function given the 
: 50% point th, the slope q, and the amplitude a.
FUNCTION SigmoidRate(v,th,a,q) {
    if (fabs(v-th) > 1e-6) {
        SigmoidRate = a * (v - th) / (1 - exp(-(v - th)/q))
    } else {
        SigmoidRate = a * q
    }
}   


: ----------------------------------------------------------------
: sign trap - trap for negative values and replace with zero
FUNCTION strap(x) {
    if (x < 0) {
        strap = 0
VERBATIM
        fprintf (stderr,"skv.mod:strap: negative state");
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
    fprintf(stderr, "skv.mod:ChkProb: argument not a probability.\n");
ENDVERBATIM
  }
}
