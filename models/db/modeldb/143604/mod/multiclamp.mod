COMMENT
  multiclamp.mod
  Generates a train of symmetrically trapezoidal current changes.
  User specifies trf (duration of rise/fall ramps), tp (duration of plateau),
  del (onset of first trapezoid), per (period, i.e. interval between trapezoid onsets), 
  and number of trapezoids.
  Ensures that period is longer than trapezoid duration.
  4/1/2012 NTC modified from TrapSyn.mod on 4/2/2012 TMM
ENDCOMMENT

NEURON {
  POINT_PROCESS MultIClamp
  ELECTRODE_CURRENT i
  RANGE trf, tp
  RANGE del, per, number
  RANGE amp, i
:  NONSPECIFIC_CURRENT i
}

UNITS {
  (mV) = (millivolt)
  (nS) = (nanosiemens)
  (nA) = (nanoamp)
}

PARAMETER {
  trf (ms) <0, 1e9> : duration of rising and falling phases
  tp  (ms) <0, 1e9> : duration of plateau
  del (ms) <0, 1e9> : latency of first transient
  per (ms) <0, 1e9> : period, i.e. interval between transient onsets
  number : how many to deliver
  amp (nA) <0, 1e9> : conductance during plateau
}

ASSIGNED {
  v (mV)
  i (nA)
  on
  tally : how many more to deliver
  m (1/ms)
  b (1)
  dur (ms)
  t0 (ms)
}

INITIAL {
  if (trf <= 0) {
    trf = 0.025 : use default time step as default rise/fall time
UNITSOFF
    printf("time rise fall must be longer than 0: ")
    printf("increased to trf = %g ms\n", trf)
UNITSON
  }
  if (tp < 0) {
    tp = 0
UNITSOFF
    printf("time plateau must not be negative: ")
    printf("changed to tp = %g ms\n", tp)
UNITSON
  }
  dur = 2*trf + tp
  if (per <= dur) {
    per = dur : allows this mechanism to be used as triangle wave
UNITSOFF
    printf("period must be longer than trapezoid duration %g: ",dur)
    printf("increased to per = %g ms\n", per)
UNITSON
  }
  on = 0
  m = 0
  b = 0
  tally = number
  if (tally > 0) {
    net_send(del, 1)
    tally = tally - 1
  }
}

BREAKPOINT {
  i = amp * (m*(t-t0) + b)
}

NET_RECEIVE (w) {
  if ((on == 0) && (flag == 1)) {
    : enter rising phase
    t0 = t
    m = 1/trf
    b = 0
    on = 1
    : prepare for plateau phase
    net_send(trf, 2)
  }
  if (flag == 2) {
    : enter plateau
    m = 0
    b = 1
    : prepare for falling phase
    net_send(tp, 3)
  }
  if (flag == 3) {
    : enter falling phase
    t0 = t
    m = -1/trf
    b = 1
    : prepare to end
    net_send(trf, 4)
  }
  if (flag == 4) {
    : end
    m = 0
    b = 0
    on = 0
    if (tally > 0) {
      : prepare to turn it on again
      net_send(per - dur, 1)
      tally = tally - 1
    }
  }
}
