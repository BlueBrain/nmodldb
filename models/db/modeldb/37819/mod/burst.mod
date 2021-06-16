: $Id: burst.mod,v 1.35 2004/05/05 20:30:04 billl Exp $

NEURON { 
  ARTIFICIAL_CELL BURST
  RANGE interval, number, numbersav
  RANGE taum, refrac, m, refractory
  RANGE taua, adap, addad, fflag
  GLOBAL debug
}

PARAMETER {
  interval      = 10 (ms) <1e-9,1e9>: time between spikes (msec)
  number        = 10 <0,1e9>          : number of spikes
  taum = 10 (ms)
  taua = 100 (ms)
  addad = 0.2
  refrac = 200 (ms)
  debug = 0
  fflag           = 1             : don't change
}

ASSIGNED {
  event (ms)
  on
  end (ms)
  m
  adap
  numbersav
  t0m(ms)
  t0a(ms)
  refractory
}

INITIAL {
  on = 0
  m = 0
  adap = 1
  refractory = 0
  t0a = 0

  index = 0
  recval()
} 

NET_RECEIVE (w) {
  if (flag == 0) { : external event
    m = m*(1-(t - t0m)/taum) : linear decay of voltage
    t0m = t
    if (m>0) { m = m + w } else { m=w } : boost the voltage
    if (t0a>0) { adap = adap*(1-(t - t0a)/taua) }
    t0a = t
    if (m>1 && refractory==0) { : threshold: start burst  
      refractory=1
      on = 1
      event = t
      numbersav=number
      if (adap>1) { number=number/adap }
      end = t + 1e-6 + interval*(number-1)
      net_send(0, 1)
    }
  }
  if (flag == 1 && on == 1) { : generate a spike
    VERBATIM
    if (debug==1) {printf("a:%g,%g,%g\n",t, adap, t0a);}
    ENDVERBATIM
    if (t0a>0) { adap = adap*(1-(t - t0a)/taua) }
    t0a = t
    if (adap+addad>=1) { adap = adap + addad } else { adap=1.0 }
    net_event(t)
    event = event + interval*adap
    net_send(event - t, 1)
    if (event>end) { 
      net_send(refrac,2) 
      number=numbersav
      on = 0 
      m = 0 
    }
  }
  if (flag == 2) { refractory = 0 }
  recval()
}

INCLUDE "ppsav.inc"
INCLUDE "pointer.inc"
