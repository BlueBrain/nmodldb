: $Id: tcifb.mod,v 1.18 2003/03/11 17:05:13 billl Exp $

COMMENT
AP occurs when v>vthresh
ENDCOMMENT

NEURON {
  POINT_PROCESS IFB
  NONSPECIFIC_CURRENT i
  GLOBAL hinf,tauh,tauh1,tauh2,vh,erev,vthresh,vreset
  RANGE gmax,gshort
}

PARAMETER {
  tauh1 = 100 (ms)
  tauh2 = 20 (ms)
  vh = -60
  erev = 120
  vthresh = -35
  vreset = -50
  gmax = 0.07 (mS/cm2)
  gshort = 0
}

ASSIGNED {
  v		(mV)		: postsynaptic voltage
  g
  i
  hinf
  tauh
}

STATE { h }

INITIAL { 
  h = 1   
  gshort = 0 
}

BREAKPOINT {
  if (v<vh) { hinf=1 tauh=tauh2
  } else {    hinf=0 tauh=tauh1  }
  SOLVE states METHOD cnexp
  if (v>vh) { : Heaviside thesholding for m
    i = gmax * h * (v-erev) + gshort*(v-vreset)
  } else {
    i=0
  }
}

DERIVATIVE states {
  h' = (hinf-h)/tauh
}

NET_RECEIVE (w) {
  if (flag==0) {
    net_event(t)
    net_send(1e-4,1) : any non zero flag will do 
    state_discontinuity(gshort, 1e5) : short circuit
  } else {
    state_discontinuity(gshort, 0)   : close the short
  }
}
