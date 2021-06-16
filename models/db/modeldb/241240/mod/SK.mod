TITLE SK channel
: SK channel following Hirschberg et al 1999 Biophys J 77:1905-1913 and Hirschberg et al 1998 J Gen Physiol 111:565-581


NEURON {
	SUFFIX 	SK
	USEION 	ca READ cai, cao, ica
	USEION	k READ ek WRITE ik
	RANGE 	gbar, g, i, inf_c, tau_c
:	GLOBAL	inf_c, tau_c
}

UNITS {
	(molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
	(pS)	=	(picosiemens)
	(um)	=	(micrometer)


	:FARADAY = 96520 (coul)
	:R = 8.3134 (joule/degC)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	gbar		= 	5 	(pS/um2)  :  Maximum conductance
	
	hill_c	=	4.64	: 2
	hill_t	=	2
	K_c		=	.00056	(mM)	: 0.56 uM
	tauA_c	=	45  	(ms)  : 50
	tau0_c	=	5	(ms)
:	shat_c	=	0.00035	(mM)
:	tauG_c	=	0.67			: left-right skew (0-1)
	diff		=	1	(1)	: diffusion factor
	cai0		=	50e-6	(mM)
	scale		=	1	(1)	: scaling for diffusion
}

ASSIGNED { 
:	celsius		(degC) : 32
	v		(mV)
	i		(mA/cm2)	
	ik		(mA/cm2)
	g		(pS/um2)
	ek		(mV)
	ica		(mA/cm2)
	cai		(mM)
	cao		(mM)
	inf_c
	tau_c		(ms)

}

STATE {
	c	: calcium dependent activation
}		

BREAKPOINT {
	SOLVE states METHOD cnexp
	g 	= 0
	g 	= gbar  * c
	i 	= g * (v - ek)*(1e-4)
	ik 	= i
}

INITIAL {
	cai=cai0
	rates( cai0)
	c = inf_c
}

DERIVATIVE states {
	rates(cai0+scale*(cai-cai0)) 
	c' = ( inf_c - c) / tau_c
}

PROCEDURE rates ( cai ( mM)) {
	inf_c = 1/(1 + (K_c/cai)^hill_c)
 	tau_c = tau0_c + tauA_c/(1 + (cai/K_c)^hill_t)
:	tau_c = tau0_c + 4*sqrt(tauG_c*(1-tauG_c))*tauA_c/(exp(tauG_c*(cai - K_c)/shat_c)+exp(-(1-tauG_c)*(cai-K_c)/shat_c))
}
	

