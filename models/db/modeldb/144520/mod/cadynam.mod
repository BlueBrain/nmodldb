TITLE Cardiac intracellular calcium accumulation, buffering and uptake
COMMENT
	modified From DiFrancesco & Noble 1985 Phil Trans R Soc Lond 307:353-398 
    modified for Neuron by FE GANNIER
	francois.gannier@univ-tours.fr (University of TOURS)
ENDCOMMENT
INCLUDE "Unit.inc"
INCLUDE "Volume.inc"
NEURON {
	SUFFIX Cadynam
	USEION ca READ ica, cai WRITE cai VALENCE 2
	USEION cr READ cri WRITE cri VALENCE 2
	USEION cu READ cui WRITE cui VALENCE 2

	RANGE Irel, Iup, Itr, m_a, m_b
}

ASSIGNED {
	v (mV)
	ica		(mA/cm2)
	Iup		(mM/ms)
	Irel	(mM/ms)
	Itr		(mM/ms)
	m_a (/ms)
	m_b (/ms)
}

PARAMETER {
	cumx 	= 5			(mM)
	Kmca	= 0.001		(mM)
	tauup = 0.025		(s)
	taurel = 0.05		(s)
	tautr = 2			(s)
	p = 1
}

STATE { : p
	cai START 5e-5	(mM)
	cri START 1		(mM)  
	cui START 2		(mM)
	m
}

LOCAL ViF, VrF, VuF
INITIAL {
	VERBATIM
		cai = _ion_cai;
		cri = _ion_cri;
		cui = _ion_cui;
	ENDVERBATIM
	m = m_a /(m_a + m_b)
	ViF = (1e-3)*Vi*F/S
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit

}

DERIVATIVE state { 
	rate(v)
	m' = m_a*(1 - m) - m_b*m
	
	Iup =  (0.001)*cai*(1 - cui/cumx)/tauup
	Itr =  (0.001)*m*(cui - cri)/tautr
	Irel = (0.001)*cri * cai^2/taurel/(cai^2 + Kmca^2)
	
	cui' = Iup*(Vi/Vup) - Itr*(Vrel/Vup)
	cri' = (Itr - Irel)
	cai' = - ica/(2*ViF) - Iup + Irel*Vrel/Vi
}

PROCEDURE rate(v(mV)) { 
	m_a = (0.001)* 0.625(/mV/s)*(v + 34)/(exp((v + 34)/4(mV)) - 1)
	m_b = (0.001)*5(/s)/ (1 + exp(-(v + 34)/4(mV)))
}
