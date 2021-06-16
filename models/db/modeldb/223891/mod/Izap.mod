COMMENT
Chirp/Swept Sine/Zap function as a pointprocess. 
Code amalgamated from the Neuron board
http://www.neuron.yale.edu/phpbb/viewtopic.php?f=8&t=897
Thanks to Ted and users kelvin and crutchley.
ENDCOMMENT

NEURON {
  POINT_PROCESS Izap
  RANGE del, dur, f0, f1, amp, i
  ELECTRODE_CURRENT i
}

UNITS {
  (nA) = (nanoamp)
  PI = (pi) (1)
}

PARAMETER {
  del (ms)
  dur (ms)
  f0 (1/s)  : frequency is in Hz
  f1 (1/s)
  amp (nA)
}

ASSIGNED {
  f (1/s)
  i (nA)
}

INITIAL {
  i = 0

  if (del<0) { del=0 }
  if (dur<0) { dur=0 }
  if (f0<=0) { f0=0 (1/s) }
  if (f1<=0) { f1=0 (1/s) }

}


BREAKPOINT {
  at_time(del)
  at_time(del + dur)

  if (t < del) {
    i=0   
  } else { 
    if (t < del+dur) {
      i = amp*sin(2*PI * (t - del) * (f0 + (f1 - f0)*(t-del)/(2*dur)) * 0.001)
    } else { 
      i = 0
    }
  }
}
