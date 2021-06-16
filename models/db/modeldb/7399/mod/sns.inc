: $Id: sns.inc,v 1.29 1997/03/25 00:05:44 billl Exp $
COMMENT
USAGE: for most receptors
*****************************************************************************
    NEURON {
      POINT_PROCESS NAME
    }

    PARAMETER {
      Cdur	= 1.08	(ms)		: transmitter duration (rising phase)
      Alpha	= 1	(/ms mM)	: forward (binding) rate
      Beta	= 0.02	(/ms)		: backward (unbinding) rate
      Erev	= -80	(mV)		: reversal potential
      Deadtime = 1	(ms)		: mimimum time between release events
      GMAX = 1		(umho)		: reference conductance
      DELAY = 0         (ms)
    }

    
    INCLUDE "sns.inc"
*****************************************************************************

USAGE: for NMDA receptor
*****************************************************************************
    NEURON{ POINT_PROCESS NMDA
      RANGE B 
    }

    PARAMETER {
      mg        = 1.    (mM)     : external magnesium concentration
      Cdur	= 1.	(ms)	 : transmitter duration (rising phase)
      Alpha	= 4.	(/ms mM) : forward (binding) rate
      Beta	= 0.0067 (/ms)	 : backward (unbinding) rate 1/150
      Erev	= 0.	(mV)	 : reversal potential
      Deadtime = 1	(ms)	 : mimimum time between release events
      GMAX     = 1      (umho)   : reference conductance
      DELAY    = 0               : axonal delay
    }

    ASSIGNED { B }

    INCLUDE "sns.inc"
    : EXTRA BREAKPOINT MUST BE BELOW THE INCLUDE
    BREAKPOINT {
      rates(v)
      g = g * B : but don't really need to readjust conductance
      i = i * B : i = g*(v - Erev)
    }

    PROCEDURE rates(v(mV)) {
      TABLE B
      DEPEND mg
      FROM -100 TO 80 WITH 180
      B = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
    }
*****************************************************************************

see end for implementation comments
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  NONSPECIFIC_CURRENT i
  RANGE R, g
  RANGE Ron, Roff, synon  : accessible for debugging
  GLOBAL GMAX, DELAY, Cdur, Alpha, Beta, Erev, Deadtime, Rinf, Rtau, q10, exptemp
}
 
INCLUDE "snsarr.inc"  : array management routines

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
  (umho) = (micromho)
  (mM) = (milli/liter)
}

PARAMETER {
  celsius
  q10 = 1.
  exptemp = 37.
}

ASSIGNED {
  v		(mV)		: postsynaptic voltage
  i 		(nA)		: current = g*(v - Erev)
  g 		(umho)		: conductance
  R				: fraction of open channels, Ron + Roff
  Ron                           : activation state while syn's turned on
  Roff                          : activation state for decaying syns
  Rinf				: steady state channels open
  Rtau		(ms)		: time constant of channel binding
  synon                         : number of syns turned on at a time
  drive                         : drive for ODE toward Rinf
  qq10                          : rate speed-up due to Q10
  edt                           : rate factor for Ron
  edb                           : decay factor for Roff
  edc                           : rate factor for increasing Rcurr
  dt
}

INITIAL {
  : initialize GLOBAL parameters, in case user changes one
  : this repeats unnecessarily for each instance
  Rinf = Alpha / (Alpha + Beta)
  qq10 = q10^((celsius-exptemp)/10.)
  Rtau = 1 / (Alpha + Beta) / qq10
  edt = exp(-dt/Rtau)
  edb = exp(-Beta * dt)
  edc = exp(-Cdur/Rtau)
  drive = Rinf * (1. - edt)

  : initialize RANGE parameters
  synon = 0
  R = 0
  Ron = 0
  Roff = 0

  : do not initialize QUEU if it hasn't been allocated by init_arrays()
  if (nsyn > 0) {
    initq()   
  } else {
    VERBATIM
    printf("WARNING nsyn==0 ");
    ENDVERBATIM
  }
}

BREAKPOINT {
    SOLVE release
    R = Ron + Roff
    g = GMAX * R
    i = g*(v - Erev)
}

PROCEDURE release() { 
  if (nsyn>0) { : do not try accessing a queue that hasn't been allocated
  VERBATIM 
  int who;
  QueU *pqueu;
  SynS *ppst;

  pqueu = (QueU *)((unsigned long) queu);

  while (t >= pqueu[(int)begsyn].time) { /*  somebody spiked delay time ago */
    ppst = (SynS *)((unsigned long) lpst);
    who = pqueu[(int)begsyn].index; /* who spiked? */
    /* calculate the decay that occurred since last activity ended */
    ppst[who].Rcurr *= exptable(-Beta*(t-ppst[who].last));
    /* transfer the value from Roff to Ron */
    Ron += ppst[who].Rcurr;
    Roff -= ppst[who].Rcurr;
    synon += ppst[who].pgm;	/*  weight syn by percent gmax */
    ppst[who].last = t + Cdur;   /* time when syn will turn off */
    popqh1(Cdur);		/* next (also add Cdur to value on the queu) */
  }

  while (t >= pqueu[(int)endsyn].time) { /*  somebody needs to be turned off */
    ppst = (SynS *)((unsigned long) lpst); /* will do this assign twice in rare cases */
    who = pqueu[(int)endsyn].index;   /* who spiked? */
    /* solve Rcurr forward in time till end of syn activity */
    ppst[who].Rcurr = ppst[who].pgm*Rinf*(1.-edc) + ppst[who].Rcurr*edc;
    Ron -= ppst[who].Rcurr;
    Roff += ppst[who].Rcurr;  /* transfer from on to off */
    synon -= ppst[who].pgm;   /* remove this percent of gmax */
    if (synon<1e-11 && synon>-1e-11) { synon=0.; }
    if (synon==0. || Ron < 0.) { Ron = 0.; }
    popqh2();  /* next */
  }

  /*  update R */
  if (synon > 0) {		/*  one or more synapses turned on? */
    Ron = synon*drive + Ron * edt; /*  drive R toward Rinf */
  } 
  Roff *= edb;			/*  Roff always decays toward 0 */
  return 0;
  ENDVERBATIM
  }
}

FUNCTION getR(x) {
  VERBATIM {
    SynS pst;
    double rnow;
    double end, rinf;

    pst = (PSTCAST[(int)_lx]);
    end = pst.last;
    rinf = pst.pgm * Rinf;

    if (end < 0.) {		/* not yet turned on */
      rnow = 0.;
    } else if (t >= end - dt/2) {	/* decay */
      rnow = pst.Rcurr * exptable(-Beta*(t-end));
    } else {			/* turning on */
      rnow = rinf + (pst.Rcurr - rinf)*exptable((t-(end-Cdur))/(-Rtau));
    }
    if (pst.pgm != 0.) {
      _lgetR = rnow/pst.pgm;  /* scale it to 1 so it looks like a state variable */
    } else {
      _lgetR = 0.;
    }
  }
  ENDVERBATIM
}

: only gets called for negative numbers
FUNCTION exptable(x) { 
  TABLE  FROM -10 TO 0 WITH 1000
  
  if ((x > -10) && (x < 0)) {
    exptable = exp(x)
  } else {
    exptable = 0.
  }
}

