TITLE Calyx of Held

COMMENT
-----------------------------------------------------------------------------

 Model of vesicle mobilization and release at multiple
 release sites in the calyx of Held.

 - basic enhanced replenishment model

 - each release site produces a pulse of transmitter, T, when
   a vesicle is released
    - T used as the pointer to the AMPA receptor

B. Graham, Dept. of Computing Science & Maths, University of Stirling
(Contact: b.graham@cs.stir.ac.uk)
(previously IANC, Division of Informatics, University of Edinburgh)

CNS 2000 Version (19/11/02)
-----------------------------------------------------------------------------
ENDCOMMENT

DEFINE SSIZE 501
DEFINE RSIZE 500

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    POINT_PROCESS COH
    RANGE T, ntot, Ttot, spike
    GLOBAL KC, KO
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (mM) = (milli/liter)
}

PARAMETER {
    dt    (ms)
    rseed = 0   : random number seed
    pv0 = 0.5  (1)  : approximate probability of release (scale factor)
    n0 = 1    (1)   : initial size of RRVP at each release site
    ke = 8 (/mM /ms)   : rate of enhanced replenishment
    kd = 0.0002  (/ms)  : background depletion of rrvp
    Camp = 0.1 (mM) : amplitude of local [Ca] transient following AP
    Cres = 0 (mM)   : amplitude of local residual [Ca]
    Cdur = 1 (ms)   : duration of local [Ca] transient
    Cnamp = 0.01 (mM)   : amplitude of distant [Ca] transient
    Cnres = 0 (mM)  : amplitude of distant residual [Ca]
    Cndur = 2 (ms)  : duration of  distant [Ca] transient
    Tamp = 1 (mM)   : amplitude of transmitter pulse
    Tdur = 1  (ms)  : duration of transmitter pulse
}

ASSIGNED {
    spike[SSIZE]    (ms)    : list of spike times
    index           : index to spike times
    tspike      (ms)    : time of last spike
    trel[RSIZE] (ms)    : time of last release
    ntot        (1) : total (or mean) RRVP size
    Ttot        (mM)    : number of releases
    km      (/ms)   : background replenishment rate
    KO[2]       (/mM /ms) : gate opening rates
    KC[2]       (/ms)     : gate closing rates
    inf[2] tau[2] fac[2]
}

STATE {
    n[RSIZE]    (1) : vesicles in RRVP
    R   (1) : probability of vesicle release
    RO[2]   (1) : gates for release
    T[RSIZE]    (mM)    : pulse of neurotransmitter
    C       (mM)    : [Ca] release transient
    Cn      (mM)    : [Ca] mobilization transient
    pv[RSIZE]   (1) : scaling for release probability
    Ccount (1)  : count of time steps for Ca transient (BPG 10-1-02)
    Cncount (1)  : count of time steps for mob Ca transient (BPG 13-1-02)
    Tcnt[RSIZE] (1)  : count of time steps for T transient (BPG 13-1-02)
}

INITIAL {
    index = 0
    FROM i = 0 TO RSIZE-1 {
      n[i] = n0
      T[i] = 0
      trel[i] = 0
      pv[i] = pv0
    }
    R = 0
    C = Cres
    Cn = Cnres
    ntot = 0
    Ttot = 0
    set_seed(rseed)
    km = n0 * kd        : background replenishment rate
    KO[0] = 150 : fast gate opening
    KO[1] = 1 : slow gate opening
    KC[0] = 30 : fast gate closing
    KC[1] = 0.1 : slow gate closing
}

BREAKPOINT {
    SOLVE release
}

PROCEDURE release() {

  if (index < SSIZE && t>=spike[index] && C < Camp) {    : presynaptic spike
    C = Camp    : calcium transient
    index=index+1 : index of next spike
    tspike = t    : time of spike
    Ccount = Cdur / dt  : number of time steps for Ca transient
    Ttot = 0
  }

  prel()    : probability of release

  ntot = 0
  FROM i = 0 TO RSIZE-1 {

    if (unirand() < dt*km) {n[i] = n[i]+1}  : background replenishment
    if (n[i] > 0 && unirand() < dt*kd) {n[i] = n[i]-1}  :depletion
    if (Cn > 0) {
      if (unirand() < dt*ke*Cn) {n[i] = n[i]+1} : extra replenishment
    }
    
    onerel(i)   : single vesicle release

    if (T[i] > 0) {
      Tcnt[i] = Tcnt[i] - 1
      if (Tcnt[i] < 0) {T[i] = 0}    : end of pulse
    }

    ntot = ntot + n[i] : total of available vesicles
  }

  if (C > Cres) {
    Ccount = Ccount - 1
    if (Ccount < 0) {
      C = Cres  : end of release transient
      Cn = Cnamp  : beginning of mobilization transient
      Cncount = Cndur / dt  : number of time steps for Ca transient
    }
  }
  
  if (Cn > Cnres) {
    Cncount = Cncount - 1 
    if (Cncount < 0) {Cn = Cnres}
  }
    
  VERBATIM
  return 0;
  ENDVERBATIM
}


PROCEDURE onerel(i) {   : one release per site per AP
  if (trel[i] < tspike) {
    : probability scaled relative to dt of 1/40 msecs (31-1-01)
    if (n[i] > 0 && unirand() < R*pv[i]*n[i]*dt*40) { : scaling (31-1-01)
      n[i] = n[i] - 1     : release
      T[i] = Tamp     : pulse of transmitter
      Tcnt[i] = Tdur / dt     : duration of pulse of transmitter
      trel[i] = t     : time of release
      Ttot = Ttot + 1     : count total releases (BPG 10-1-00)
    }
  }
}


PROCEDURE prel() {  : Probability of release
  rates(C)
  R = 1
  FROM i=0 TO 1 {
    RO[i] = RO[i] + fac[i]*(inf[i] - RO[i])
    R = R * RO[i]
  }
}

 
PROCEDURE rates(C) {LOCAL a, b  :Computes gate rates at concentration C.
        TABLE inf, fac, tau DEPEND dt FROM 0 TO 0.2 WITH 200
        FROM j=0 TO 1 {
            a = KO[j] * C
            b = KC[j]
            tau[j] = 1/(a + b)
            inf[j] = a/(a + b)
            fac[j] = (1 - exp(-dt/tau[j]))
        }
}


FUNCTION unirand() {    : uniform random numbers between 0 and 1
        return(scop_random())
}
