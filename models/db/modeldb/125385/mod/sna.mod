COMMENT
----------------------------------------------------------------
This a stochastic version of the na3h5.mod of Z. Mainen in
Mainen & Sejnowski 95 (Zach Mainen, Salk Institute, 1994, zach@salk.edu).

It is a sodium channel, Hodgkin-Huxley like kinetics, based on the gates
model, assuming stochastic opening and closing.

The rates function are adapted directly from the na3h5.mod
by Zach Mainen, available through:
http://www.cnl.salk.edu/~zach/initdemo.html

The stochastic model is similar to that used by Schneidman, Freedman, 
Segev, 1998 and is as follows (ascii-art):

           = 3alpha_m =>         = 2alpha_m =>        = alpha_m =>  
    [m0h1]               [m1h1]               [m2h1]             <[m3h1]>
          <= beta_m =           <= 2beta_m =         <= 3beta_m =  
        ^                  ^                    ^                 ^  
beta_h| |          beta_h| |            beta_h| |         beta_h| |
      | |alpha_h         | |alpha_h           | | alpha_h       | | alphah  
      | |                | |                  | |               | |
      V    = 3alpha_m => V      = 2alpha_m => V    = alpha_m => V  
    [m0h0]                [m1h0]               [m2h0]              [m3h0]
           <= beta_m =          <= 2beta_m =          <= 3beta_m =  


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
    SUFFIX sna
    USEION na READ ena WRITE ina     
    GLOBAL minf, hinf, mtau, htau,am,bm,ah,bh     
    RANGE N, reff, eta, gamma, deterministic, gna
    GLOBAL P_am, P_bm, P_ah, P_bh, hinf_curve
    GLOBAL tha, thi1, thi2, qa, qi, qinf, thinf, vshift
    GLOBAL Ra, Rb, Rd, Rg
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

    thi1  = -50 (mV)        : v 1/2 for inact   
    thi2  = -75 (mV)        : v 1/2 for inact   
    qi   = 5    (mV)            : inact tau slope
    hinf_curve = 1          : 0 if use normal relation for hinf
                        : otherwise use curve specified by following 
                        : two parameters
    thinf  = -65    (mV)        : inact inf slope   
    qinf  = 6.2 (mV)        : inact inf slope
    Rg   = 0.0091   (/ms)   : inact (v) 
    Rd   = 0.024    (/ms)   : inact recov (v) 

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
     ah     (/ms)
    bh      (/ms)
    minf        
     hinf
    mtau (ms)    
     htau (ms)
     tadj

     N 
     reff   (pS/um2)
    scale_dens (pS/um2)
    P_am        : probability of one channel making am transition
    P_bm
    P_ah
    P_bh
    
    wflag
}
 
STATE { m h                             : deterministic variables
        m0h0 m0h1 m1h0 m1h1 m2h0 m2h1 m3h0 m3h1 : states
m0h0_m1h0  m1h0_m2h0  m2h0_m3h0  m0h1_m1h1  m1h1_m2h1  m2h1_m3h1  
m3h0_m2h0  m2h0_m1h0  m1h0_m0h0  m3h1_m2h1  m2h1_m1h1  m1h1_m0h1  
m0h0_m0h1 m0h1_m0h0 m1h0_m1h1 m1h1_m1h0 m2h0_m2h1 m2h1_m2h0 m3h0_m3h1 m3h1_m3h0 
}
: ----------------------------------------------------------------
: Initialization.
INITIAL { 
    trates(v+vshift) 
    wflag = 1   : only give a warning once!   
    m = minf
    h = hinf
    scale_dens = gamma/area
    reff = gamma*eta
    N   = floor(eta*area + 0.5)   : round to nearest number of channels

    m1h0 = floor(3*m*(1-m)*(1-m)*(1-h)*N + 0.5)
    m2h0 = floor(3*m*m*(1-m)*(1-h)*N + 0.5)
    m3h0 = floor(m*m*m*(1-h)*N + 0.5)

    m0h1 = floor((1-m)*(1-m)*(1-m)*h*N + 0.5)
    m1h1 = floor(3*m*(1-m)*(1-m)*h*N + 0.5)
    m2h1 = floor(3*m*m*(1-m)*h*N + 0.5)
    m3h1 = floor(m*m*m*h*N + 0.5)
    
    : put the rest of the channels in the non-conducting & inactivated state
    m0h0 = N - (m1h0 + m2h0 + m3h0 + m0h1 + m1h1 + m2h1 + m3h1)

    m0h0_m1h0=0 
    m1h0_m2h0=0 
    m2h0_m3h0=0 
    m0h1_m1h1=0
    m1h1_m2h1=0
    m2h1_m3h1=0 
    m3h0_m2h0=0 
    m2h0_m1h0=0
    m1h0_m0h0=0 
    m3h1_m2h1=0 
    m2h1_m1h1=0 
    m1h1_m0h1=0

    m0h0_m0h1=0 
    m0h1_m0h0=0 
    m1h0_m1h1=0 
    m1h1_m1h0=0 
    m2h0_m2h1=0 
    m2h1_m2h0=0 
    m3h0_m3h1=0 
    m3h1_m3h0=0
}

: ----------------------------------------------------------------
: Breakpoint for each integration step
BREAKPOINT {
    SOLVE states
    if (deterministic) { 
        if (deterministic-1){        
    gna = m*m*m*h*reff*tadj    
    } else {                            
    gna = floor(m*m*m*h* N + 0.5) * scale_dens *tadj}
    } else{                                        
    gna = strap(m3h1) * scale_dens * tadj
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
    h = h + (1 - exp(-dt/htau)) * (hinf-h)
    
    P_am = strap(am*dt)
    P_bm  = strap(bm*dt)
    
    : check that will represent probabilities when used
    ChkProb( 3.0 * P_am)
    ChkProb( 3.0 * P_bm)
    ChkProb( P_bm/(1.0-2.0*P_am) )
    ChkProb( 2.0 * P_bm/(1.0-P_am) )
    
    : m gate transitions
:    if (deterministic) {
:    m0h0_m1h0 = 3.0*P_am*m0h0
:    m1h0_m2h0 = 2.0*P_am*m1h0
:    m1h0_m0h0 = P_bm/(1.0-2.0*P_am)*(m1h0 - m1h0_m2h0)  
:    m2h0_m3h0 = P_am*m2h0
:    m2h0_m1h0 = 2.0*P_bm/(1.0-P_am)*(m2h0 - m2h0_m3h0)
:    m3h0_m2h0 = 3.0*P_bm*m3h0
:    m0h1_m1h1 = 3.0*P_am*m0h1    
:    m1h1_m2h1 = 2.0*P_am*m1h1
:    m1h1_m0h1 = P_bm/(1.0-2.0*P_am)*(m1h1 - m1h1_m2h1)
:    m2h1_m3h1 = P_am*m2h1
:    m2h1_m1h1 = 2.0*P_bm/(1.0-P_am)*(m2h1 - m2h1_m3h1)
:    m3h1_m2h1 = 3.0*P_bm*m3h1
:    }
:    else {
    m0h0_m1h0 = BnlDev_RNG(3.0*P_am,m0h0)
    m1h0_m2h0 = BnlDev_RNG(2.0*P_am,m1h0)
    m1h0_m0h0 = BnlDev_RNG(P_bm/(1.0-2.0*P_am), m1h0 - m1h0_m2h0)  
    m2h0_m3h0 = BnlDev_RNG(P_am,m2h0)
    m2h0_m1h0 = BnlDev_RNG(2.0*P_bm/(1.0-P_am), m2h0 - m2h0_m3h0)
    m3h0_m2h0 = BnlDev_RNG(3.0*P_bm,m3h0)
    m0h1_m1h1 = BnlDev_RNG(3.0*P_am, m0h1)
    m1h1_m2h1 = BnlDev_RNG(2.0*P_am, m1h1)
    m1h1_m0h1 = BnlDev_RNG(P_bm/(1.0-2.0*P_am), m1h1 - m1h1_m2h1)
    m2h1_m3h1 = BnlDev_RNG(P_am,m2h1)
    m2h1_m1h1 = BnlDev_RNG(2.0*P_bm/(1.0-P_am), m2h1 - m2h1_m3h1)
    m3h1_m2h1  = BnlDev_RNG(3.0*P_bm,m3h1)
:    }
    : new numbers in each state after the m gate transitions
    m0h0 = m0h0 - m0h0_m1h0 + m1h0_m0h0
    m1h0 = m1h0 - m1h0_m2h0 - m1h0_m0h0  + m2h0_m1h0 + m0h0_m1h0
    m2h0 = m2h0 - m2h0_m3h0 - m2h0_m1h0  + m3h0_m2h0 + m1h0_m2h0
    m3h0 = m3h0 - m3h0_m2h0 + m2h0_m3h0

    m0h1 = m0h1 - m0h1_m1h1 + m1h1_m0h1
    m1h1 = m1h1 - m1h1_m2h1 - m1h1_m0h1 + m2h1_m1h1 + m0h1_m1h1
    m2h1 = m2h1 - m2h1_m3h1 - m2h1_m1h1 + m3h1_m2h1 + m1h1_m2h1
    m3h1 = m3h1 - m3h1_m2h1 + m2h1_m3h1

    : probabilities of making h gate transitions
    P_ah = strap(ah*dt)
    P_bh = strap(bh*dt)
    
    ChkProb(P_ah)
    ChkProb(P_bh)
    
    : number making h gate transitions
:    if (deterministic) {
:    m0h0_m0h1 = P_ah*m0h0
:    m0h1_m0h0 = P_bh*m0h1
:    m1h0_m1h1 = P_ah*m1h0
:    m1h1_m1h0 = P_bh*m1h1
:    m2h0_m2h1 = P_ah*m2h0
:    m2h1_m2h0 = P_bh*m2h1
:    m3h0_m3h1 = P_ah*m3h0
:    m3h1_m3h0 = P_bh*m3h1
:    }
:    else {
    m0h0_m0h1 = BnlDev_RNG(P_ah,m0h0)
    m0h1_m0h0 = BnlDev_RNG(P_bh,m0h1)
    m1h0_m1h1 = BnlDev_RNG(P_ah,m1h0)
    m1h1_m1h0 = BnlDev_RNG(P_bh,m1h1)
    m2h0_m2h1 = BnlDev_RNG(P_ah,m2h0)
    m2h1_m2h0 = BnlDev_RNG(P_bh,m2h1)
    m3h0_m3h1 = BnlDev_RNG(P_ah,m3h0)
    m3h1_m3h0 = BnlDev_RNG(P_bh,m3h1)
:    }
    m0h0 = m0h0 - m0h0_m0h1  + m0h1_m0h0
    m1h0 = m1h0 - m1h0_m1h1  + m1h1_m1h0
    m2h0 = m2h0 - m2h0_m2h1  + m2h1_m2h0
    m3h0 = m3h0 - m3h0_m3h1  + m3h1_m3h0

    m0h1 = m0h1 - m0h1_m0h0  + m0h0_m0h1
    m1h1 = m1h1 - m1h1_m1h0  + m1h0_m1h1
    m2h1 = m2h1 - m2h1_m2h0  + m2h0_m2h1
    m3h1 = m3h1 - m3h1_m3h0  + m3h0_m3h1
}


: ----------------------------------------------------------------
: trates - compute rates, using table if possible
PROCEDURE trates(vm) {     TABLE minf, mtau, hinf, htau, am, bm, ah, bh, tadj
    DEPEND dt,Ra,Rb,Rd,Rg,tha,thi1,thi2,qa,qi,qinf,q10,orig_temp,celsius, hinf_curve
    FROM vmin TO vmax WITH 199
        tadj = q10^((celsius-orig_temp)/10)
    
    : m activation variable
    am = SigmoidRate(vm,tha,Ra,qa)
    am = am * tadj
    bm = SigmoidRate(-vm,-tha,Rb,qa)
    bm = bm * tadj
    mtau = 1/(am+bm)
    minf = am*mtau
    
    : h inactivation variable
    ah = SigmoidRate(vm,thi1,Rd,qi)
    ah = ah * tadj
    bh = SigmoidRate(-vm,-thi2,Rg,qi)
    bh = bh * tadj
    htau = 1/(ah+bh)
    : hinf_curve is non-zero if using explicit fit, zero otherwise
    if (hinf_curve == 0) {
        hinf = ah*htau
    }
    else {
        hinf = 1/(1+exp((vm-thinf)/qinf))
        : recompute rates based on hinf
        ah = hinf/htau
        bh = 1/htau - ah
    }
    
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
        fprintf (stderr,"sna.mod:strap: negative value for state");
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
    fprintf(stderr, "sna.mod:ChkProb: argument not a probability.\n");
ENDVERBATIM
    wflag =0}
  } }
