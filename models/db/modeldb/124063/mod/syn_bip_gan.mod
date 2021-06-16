: Graded Synapse with first order binding kinetics

NEURON {
POINT_PROCESS GradSyn_bip_gan
POINTER V_pre
RANGE e, tau, s_inf, V_thr, V_slope, g_max
NONSPECIFIC_CURRENT i
}

PARAMETER {
	tau = 10(ms)
	e = -70 (millivolts)
	g_max = 0.00256 (umho)	: maximal conductance dAMPA
	V_slope = 10
	V_thr = -38.5
}

ASSIGNED {
	v (millivolt)
	V_pre (millivolt)
	i (nanoamp)
	s_inf
	g
}

STATE {
	s 
}
 
BREAKPOINT {
	SOLVE state METHOD euler
	g = g_max * s
	i = g * (v - e)
}
 
INITIAL {
	s = 0.1
}
 
  
DERIVATIVE state {
	if (V_pre <= V_thr) {
		s_inf = 0
	}
	
	if (V_pre > V_thr){
		s_inf = tanh((V_pre-V_thr)/V_slope)
	}
	
	if (s < 1e-06) {
		s = 1e-06
	} 
	
	if (s > 1) {
		s = 1
	}
	
	s'=(s_inf-s)/((1-s_inf)*tau*s)
}
