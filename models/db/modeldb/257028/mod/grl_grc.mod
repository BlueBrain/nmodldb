TITLE Granule cell in the granular layer (GrL)
COMMENT
    Modified from Garrido et al, 2013.
ENDCOMMENT 

NEURON {
	POINT_PROCESS GrC
	RANGE v_GrC, g_ampa, tau_ampa, e_ampa, i_ampa, g_nmda, tau_nmda, e_nmda, i_nmda, g_gaba, tau_gaba, e_gaba, i_gaba
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	Cm = 2 (pF)
	epas = -65 (mV)
	Grest = 0.2 (nS)
	
	e_ampa = 0	(mV)
	tau_ampa = 0.5 (ms)
	
	e_nmda = 0	(mV)
	tau_nmda = 40 (ms)
	alp = 62e-3 (1/mV)
    bet = 3.57 (mM)
	MgC = 1.2 (mM)
	
	e_gaba = -65 (mV)
	tau_gaba = 10 (ms)
}

ASSIGNED {
	i_ampa (nA)
	i_nmda (nA)
	i_gaba (nA)
}

STATE {
	v_GrC (mV)
	g_ampa (uS)
	g_nmda (uS)
	g_gaba (uS)
}

INITIAL {
	v_GrC = epas
	g_ampa=0
	g_nmda=0
	g_gaba=0
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	if (v_GrC>10 && v_GrC<20) {
		v_GrC = epas
		}
	i_ampa = g_ampa*(e_ampa - v_GrC)
	i_nmda = g_nmda*(e_nmda - v_GrC)/(1+exp(-alp*v_GrC*MgC/bet))
	i_gaba = g_gaba*(e_gaba - v_GrC)
}

DERIVATIVE state {
	v_GrC' = (i_ampa + i_nmda + i_gaba + Grest*(epas - v_GrC))/Cm
	g_ampa' = -g_ampa/tau_ampa
	g_nmda' = -g_nmda/tau_nmda
	g_gaba' = -g_gaba/tau_gaba
}

NET_RECEIVE(weight (uS)) {
	if (weight>=1 && weight<=5) {
		g_ampa = g_ampa + weight
	}
	if (weight>=0.1 && weight<=0.8) {
		g_nmda = g_nmda + 0.087*4
	}
	if (weight>=5) {
		g_gaba = g_gaba + weight
	}
	: Spike detection; spike if membrane potential>-40 mV and given enough input
	if (weight>=0.01 && weight<=0.03 && v_GrC>-40) {
		v_GrC = 20
		net_event(t) : Release a spike
		g_ampa=0
		g_nmda=0
		g_gaba=0
	}
}
