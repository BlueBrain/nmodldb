TITLE is :Inhibitory synapse in the pyramidal cell

NEURON {	
	POINT_PROCESS is
        RANGE tau1, tau2, i, g
	USEION cl READ ecl WRITE icl VALENCE -1	
        		   
	POINTER e
}

PARAMETER {
	tau1 = 2   (ms)
	tau2 = 6   (ms)
}

ASSIGNED {
	v    (mV)
	e    (mV)
	ecl  (mV)
	icl  (nanoamp)
  g (microsiemens)
  factor
}

STATE {
   A (microsiemens)
   B (microsiemens)
}

INITIAL {
 LOCAL tp
 if (tau1/tau2 > .9999) {
     tau1 = .9999*tau2
 }
 A = 0
 B = 0
 tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
 factor = -exp(-tp/tau1) + exp(-tp/tau2)
 factor = 1/factor
}

BREAKPOINT {
 SOLVE state METHOD cnexp
 g = B - A
 icl = g*(v - e)
}

DERIVATIVE state {
 A' = -A/tau1
 B' = -B/tau2
}

 NET_RECEIVE(weight (microsiemens)) {
 A = A + weight*factor
 B = B + weight*factor
} 
