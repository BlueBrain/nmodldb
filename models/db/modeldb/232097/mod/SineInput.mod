
NEURON {
  POINT_PROCESS IClampSine
  RANGE i,del,dur,u0,u1,r,torn,std,bias,A,f
  ELECTRODE_CURRENT i
}

UNITS {
  (nA) = (nanoamp)
}

PARAMETER {
  del = 100    (ms)
  dur = 200   (ms)
:  torn= 500   (ms)
:  std = 0.2   (nA)
:  u0  = 0.2   (nA)
:  u1  = 0.8   (nA)
:  r   = 60 
  bias = 0    (nA)
  pi = 3.1415926
  A = 1 (nA)
  f = 2.5  : frequency of sine
}

ASSIGNED {
  ival (nA)
  i (nA)
  amp (nA)
  noise (nA)
  res   (nA) :Respiration modulation
  on (1)
}

INITIAL {
  i  = 0
  on = 0
  net_send(del, 1)
}

BEFORE BREAKPOINT {
  if  (on) {
    :noise = normrand(0,std*1(/nA))*1(nA)
	res = 0.5*sin(2*pi*f*(t/1000-0.25/f))+0.5
    :amp = u0 + 0.5*(u1-u0)*(tanh((t-torn)/(r/3)/(1(ms))-3)+1)
    
	ival = A*res + bias
  } else {
    ival = 0
  }
}

BREAKPOINT {
  i = ival
}

NET_RECEIVE (w) {
  if (flag == 1) {
    if (on == 0) {
      : turn it on
      on = 1
      : prepare to turn it off
      net_send(dur, 1)
    } else {
      : turn it off
      on = 0
    }
  }
}


