NEURON {
	SUFFIX ab4BK2
	USEION ca READ cai
	USEION k READ ek WRITE ik
	RANGE gbar,  g, i
	GLOBAL	inf_n, tau_n
}
UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(mA) = (milliamp)
	(pS)	= (picosiemens)
	(um)	= (micrometer)
	(S) = (siemens)
:	B= 0.26 (mM-cm2/mA-ms)
}
PARAMETER {
	gbar 	= 3	(pS/um2)
	
	:  inf_n parameters
	vhalf_n 	= -56.5 	(mV)
	slope_n 	= 11.8	(mV)

	:  tau_n parameters
	tauA_n	= 14.5		(ms)
	tauG_n	= 0.48	: Left-right bias. range(0,1)
	tau0_n	= 1		(ms)  : minimum tau

	:  calcium influence on V0.5
	bv 	= 	0.003		(mM)
	nv 	= 2
	ah 	= 120		(mV)

	:  calcium influence on tau
	bt	=	0.0022		(mM)
	nt	=	4
	th	=	140		(ms)

	:  calcium conc parameters
	scale	= 1	(1)	:scale cai for nanodomain
	cai0 	= .00005	(mM)

}
ASSIGNED {
	v		(mV)
	celsius	(degC)
	ek		(mV)
	ik		(mA/cm2)
	ica		(mA/cm2)
	i		(mA/cm2)
	area		(microm2)
	g		(pS/um2)
	inf_n		
	tau_n		(ms)
	vh		(mV)
	ta		(ms)
 	cai		(mM)
:	cal 		(mM)
}
:STATE { ca_i (mM)  n }
STATE { n }

BREAKPOINT {
	SOLVE state METHOD cnexp
	g=0
	if( n >=0 ) {  
	   g = gbar * n
	}
	i = g * (v - ek) * (1e-4)
	ik = i
 
}

INITIAL {
	rates(v, cai)
	n = inf_n
}

DERIVATIVE state {
:	ca_i' = -B*(ica*cascale) - (ca_i - ca0)/tau
	rates (v, cai0+scale*(cai-cai0))
	n' = (inf_n - n) / tau_n
}


PROCEDURE rates (v (mV), cai (mM)) {

	vh 	= vhalf_n + ah /( 1+(cai/bv)^nv ) : Ca shift in v0.5 

	inf_n	= 1 / (1 + exp((v-vh)/-slope_n))

	ta	= tauA_n + th / ( 1 + (cai/bt)^nt) : Ca shift in tauA

	tau_n	= tau0_n + 4*ta*sqrt(tauG_n*(1-tauG_n))/(exp(tauG_n*(v - vh)/slope_n)+exp(-(1-tauG_n)*(v - vh)/slope_n))

}
