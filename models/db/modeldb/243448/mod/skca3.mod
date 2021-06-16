TITLE SKCa current for bladder small DRG neuron soma model
: Author: Darshan Mandge (darshanmandge@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: A biophysically detailed computational model of bladder small DRG neuron soma 
: Darshan Mandge and Rohit Manchanda, PLOS Computational Biology (2018)

UNITS {
	(molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
	FARADAY = (faraday)  (kilocoulombs)
	R = (k-mole) (joule/degC)
	(nA)=(nanoamp)
	(um)=(micrometer)
}

NEURON {
	SUFFIX skca3
	USEION ca READ cai
	USEION k READ ek WRITE ik
	RANGE ik, gbar, g
	RANGE oinf, hcsk3, E50hsk3
	RANGE m_vh, m_sf, m
	THREADSAFE
}

PARAMETER {
	v				(mV)
	gbar = 0.0009	(mho/cm2)
	cai				(mM) 
	ek				(mV)
		
	: Strobeak 2006
	hcsk3	= 5.6 	(1)
	E50hsk3 = 0.42e-3 (mM)

	m				(1)
	m_vh = 24		(mV) 
	m_sf = 128		(mV)
}

ASSIGNED {
	ik		(mA/cm2)
	o
	oinf
	g		(mho/cm2)
}

BREAKPOINT {
	rate(cai,v)
	g = gbar*o*m
	ik = gbar*o*m*(v - ek)
}

INITIAL {
	rate(cai,v)
}

FUNCTION_TABLE tabvh(cai(mM)) (mV) 
FUNCTION_TABLE tabsf(cai(mM)) (mV) 

PROCEDURE rate(ca (mM), v(mV)) {
	UNITSOFF
	: 1. caclium dependent activation (No time dependence)
	: Strobaek et al., 2006
	oinf = (ca^hcsk3)/((E50hsk3^hcsk3)+(ca^hcsk3))
	UNITSON
	o=oinf
	
	: 2. Voltage dependent activation (rectification)
	: Data Fits: Hougaard et al., 2009 and Strobaek et al., 2006
	m_vh = tabvh(ca)
	m_sf = tabsf(ca)
	
	m = 1/(1+exp((v-(ek+m_vh))/(m_sf))): voltage dependent activation
}