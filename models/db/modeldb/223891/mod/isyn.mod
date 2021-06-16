NEURON {
	POINT_PROCESS isyn
	RANGE del, amp, tau1, tau2, i, factor
	ELECTRODE_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
}

PARAMETER {
	del=0 (ms)
	tau1=.5 (ms)	<1e-3,1e6>
	tau2=1 (ms)   <1e-3,1e6>
	amp=0 	(nA)	<0,1e9>
	factor
}

ASSIGNED {
  v (mV)
  i (nA) 
}

INITIAL {
  if (tau1==tau2) {
    tau1 = tau1*0.999
  }
  factor = -1*((tau2/tau1)^(tau2/(tau1-tau2)))*((tau2-tau1)/tau1)
}

BREAKPOINT {
	if (amp) { at_time(del) }
	i = twoexp( (t - del) )
}

FUNCTION twoexp(x(ms)) (nA) {
  if (x < 0 || x/tau2 > 10) {
    twoexp = 0
  }else{
    twoexp = -amp/factor*(exp(-x/tau2)-exp(-x/tau1))
  }
}

