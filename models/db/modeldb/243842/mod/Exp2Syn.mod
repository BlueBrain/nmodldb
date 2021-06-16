: based on AlphaSynapse and Exp2Syn
: like AlphaSynapse, produces a single conductance change at a predetermined time
: like Exp2Syn, time course of conductance change is difference between two exponentials
: tau1 must be < tau2
: if user specifies tau1 >= tau2, initialization will force tau1 to be slightly < tau2,
: and time course will closely approximate an alpha synapse with tau = tau2
: Matt McGinley and Ted Carnevale, 7/21/2010

NEURON {
   POINT_PROCESS Exp2Delay
   RANGE onset, tau1, tau2, gmax, e, i
   NONSPECIFIC_CURRENT i
   RANGE g
}

UNITS {
   (nA) = (nanoamp)
   (mV) = (millivolt)
   (uS) = (microsiemens)
}

PARAMETER {
   onset = 0 (ms)
   tau1 = 0.1 (ms) <1e-9,1e9>
   tau2 = 1.5 (ms) <1e-9,1e9>
   gmax = 0 (uS)
   e = 0   (mV)
}

ASSIGNED { v (mV) i (nA) g (uS) factor }

INITIAL {
   LOCAL tp
   if (tau1/tau2 > .9999) {
      tau1 = .9999*tau2
   }
   tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
   factor = -exp(-tp/tau1) + exp(-tp/tau2)
   factor = 1/factor
}

BREAKPOINT {
   if (gmax) { at_time(onset) }
        g = gmax*dblexp(t-onset)

   i = g*(v - e)
}

FUNCTION dblexp(x (ms)) {
  if (x < 0 || x/tau2 > 10) {
    dblexp = 0
  } else {
    dblexp =  0.2 *factor*(-exp(-x/tau1) + exp(-x/tau2))
  }
}