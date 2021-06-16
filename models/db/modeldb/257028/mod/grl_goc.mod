TITLE Golgi cell in the granular layer (GrL)
COMMENT
    Modified from Garrido et al, 2013.
ENDCOMMENT 

NEURON {
	POINT_PROCESS GoC
	RANGE v_GoC, g_ampa, tau_ampa, e_ampa, i_ampa, g_gaba, tau_gaba, e_gaba, i_gaba
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	Cm = 50 (pF)
	epas = -65 (mV)
	Grest = 3 (nS)
	
	e_ampa = 0	(mV)
	tau_ampa = 0.5 (ms)
	
	e_gaba = -65 (mV)
	tau_gaba = 10 (ms)
}

ASSIGNED {
	i_ampa (nA)
	i_gaba (nA)
}

STATE {
	v_GoC (mV)
	g_ampa (uS)
	g_gaba (uS)
}

INITIAL {
	v_GoC = epas
	g_ampa=0
	g_gaba=0
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	if (v_GoC>10 && v_GoC<20) {
		v_GoC = epas
		}
	i_ampa = g_ampa*(e_ampa - v_GoC)
	i_gaba = g_gaba*(e_gaba - v_GoC)
}

DERIVATIVE state {
	v_GoC' = (i_ampa + i_gaba + Grest*(epas - v_GoC))/Cm
	g_ampa' = -g_ampa/tau_ampa
	g_gaba' = -g_gaba/tau_gaba
}

NET_RECEIVE(weight (uS)) {
	if (weight>=15) {
		g_ampa = g_ampa + weight
	}
	if (weight>=1 && weight<=15) {
		g_gaba = g_gaba + weight
	}
	: Spike detection; spike if membrane potential>-50 mV and given enough input
	if (weight>=0.01 && weight<=0.03 && v_GoC>-50) {
		v_GoC = 20
		net_event(t) : Release a spike
		g_ampa=0
		g_gaba=0
	}
}
