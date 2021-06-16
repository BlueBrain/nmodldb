COMMENT
----------------------------------------------------------------

Persistent version of sna.mod (inactivation removed)

The model keeps track of the number of channels in each state and 
uses a binomial rng to update this number. The model requires that the 
RNG mechanism be inserted somewhere in the simulation in order to provide
the BnlDev_RNG function.

Jan 1999 Mickey London, Hebrew University, 1998, mikilon@lobster.ls.huji.ac.il
        Peter N. Steinmetz, Caltech, peter@klab.caltech.edu
14 Sep 99 PNS. Added deterministic flag.   
01 Sep 02 K. Diba changed deterministic to RANGE variable     
----------------------------------------------------------------
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX snap
    USEION na READ ena WRITE ina     
    GLOBAL minf, mtau,am,bm     
    RANGE N, reff, eta, gamma, deterministic, gna
    GLOBAL P_am, P_bm
    GLOBAL tha, qa, vshift
    GLOBAL Ra, Rb
    GLOBAL vmin, vmax, q10, orig_temp, wflag, tadj
    GLOBAL DONT_VECTORIZE   : prevent vectorization to agree with RNG.mod
}
  
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (um) = (micron)
} 

PARAMETER {
    v       (mV)
    dt      (ms)
    area

    vshift = -10  (mV)        : voltage shift (affects all)
                                
    tha  = -35  (mV)        : v 1/2 for act     (-42)
    qa   = 9    (mV)            : act slope     
    Ra   = 0.182    (/ms)   : open (v)      
    Rb   = 0.124    (/ms)   : close (v)     

    gamma  =  20      (pS)    
    eta     = 50      (1/um2) : channel density
    
    celsius (degC)
    orig_temp = 23 (degC)   : original temperature
    q10 = 2.3               : temperature sensitivity
    
    deterministic = 0   : if non-zero, use deterministic variables
    
    vmin = -120 (mV)        : range to construct tables for
    vmax = 100  (mV)
    DONT_VECTORIZE          : required declaration
}

ASSIGNED {
    ina         (mA/cm2)
    gna     (pS/um2)
    ena     (mV)

     am     (/ms)
    bm      (/ms)
    minf 
    mtau (ms)
     tadj

     N 
     reff   (pS/um2)
    scale_dens (pS/um2)
    P_am        : probability of one channel making am transition
    P_bm
    
    wflag
}
 
STATE { m                             : deterministic variables
    m00 m1 m2 m3 : states
    m00_m1  m1_m2  m2_m3 
    m3_m2  m2_m1  m1_m00 
}
: ----------------------------------------------------------------
: Initialization.   

INITIAL { 
    trates(v+vshift) 
    wflag = 1   : only give a warning once!   
    m = minf
    scale_dens = gamma/area
    reff = gamma*eta
    N   = floor(eta*area + 0.5)   : round to nearest number of channels

    m1 = floor(3*m*(1-m)*(1-m)*N + 0.5)
    m2 = floor(3*m*m*(1-m)*N + 0.5)
    m3 = floor(m*m*m*N + 0.5)

    : put the rest of the channels in the non-conducting & inactivated state
    m00 = N - (m1 + m2 + m3)


    m00_m1 = 0  
    m1_m2 = 0
    m2_m3 = 0
    m3_m2 = 0 
    m2_m1 = 0 
    m1_m00 = 0  
}

: ----------------------------------------------------------------
: Breakpoint for each integration step
BREAKPOINT {
    SOLVE states
    if (deterministic) { 
        if (deterministic-1){        
    gna = m*m*m*reff*tadj    
    } else {                            
    gna = floor(m*m*m*N + 0.5) * scale_dens *tadj}
    } else{                                        
    gna = strap(m3) * scale_dens * tadj
    }
    ina = (1e-4) * gna * (v - ena)
} 


: ----------------------------------------------------------------
: states - compute state variables
PROCEDURE states() {

VERBATIM
    extern double BnlDev_RNG();
ENDVERBATIM     
    trates(v+vshift)

    : deterministic versions of state variables
    : integrated by relaxing toward the steady state value
    m = m + (1 - exp(-dt/mtau)) * (minf-m)
    
    P_am = strap(am*dt)
    P_bm  = strap(bm*dt)
    
    : check that will represent probabilities when used
    ChkProb( 3.0 * P_am)
    ChkProb( 3.0 * P_bm)
    ChkProb( P_bm/(1.0-2.0*P_am) )
    ChkProb( 2.0 * P_bm/(1.0-P_am) )

    m00_m1 = BnlDev_RNG(3.0*P_am,m00)
    m1_m2 = BnlDev_RNG(2.0*P_am,m1)
    m1_m00 = BnlDev_RNG(P_bm/(1.0-2.0*P_am), m1 - m1_m2)  
    m2_m3 = BnlDev_RNG(P_am,m2)
    m2_m1 = BnlDev_RNG(2.0*P_bm/(1.0-P_am), m2 - m2_m3)
    m3_m2 = BnlDev_RNG(3.0*P_bm,m3)
    
    : new numbers in each state after the m gate transitions
    m00 = m00 - m00_m1 + m1_m00
    m1 = m1 - m1_m2 - m1_m00  + m2_m1 + m00_m1
    m2 = m2 - m2_m3 - m2_m1  + m3_m2 + m1_m2
    m3 = m3 - m3_m2 + m2_m3
}


: ----------------------------------------------------------------
: trates - compute rates, using table if possible
PROCEDURE trates(vm) {     
    TABLE minf, mtau, am, bm, tadj
    DEPEND dt,Ra,Rb,tha,qa,q10,orig_temp,celsius
    FROM vmin TO vmax WITH 199
        tadj = q10^((celsius-orig_temp)/10)
    
    : m activation variable
    am = SigmoidRate(vm,tha,Ra,qa)
    am = am * tadj
    bm = SigmoidRate(-vm,-tha,Rb,qa)
    bm = bm * tadj
    mtau = 1/(am+bm)
    minf = am*mtau
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
: sign trap - trap negative numbers and replace with zero
FUNCTION strap(x) {
    if (x < 0) {
        strap = 0
    if (wflag){            
VERBATIM
        fprintf (stderr,"snap.mod:strap: negative value for state");
ENDVERBATIM
    wflag = 0}
    } else { 
        strap = x
    }
}

: ----------------------------------------------------------------
: ChkProb - Check that number represents a probability
PROCEDURE ChkProb(p) {
  if (p < 0.0 || p > 1.0) {
    if (wflag){
VERBATIM
    fprintf(stderr, "snap.mod:ChkProb: argument not a probability.\n");
ENDVERBATIM
    wflag =0}
  } }
