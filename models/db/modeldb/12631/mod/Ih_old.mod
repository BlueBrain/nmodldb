: $Id: Ih_old.mod,v 1.6 1995/02/16 22:18:58 ethomas Exp $
TITLE anomalous rectifier channel
COMMENT
:
: Anomalous Rectifier Ih - cation (Na/K) channel
: Differential equations
:
: Model of double activation (Destexhe & Babloyantz, 1992)
: Activation functions were fitted from 
: McCormick & Pape,  J. Physiol. 431: 291, 1990.
: and Soltesz et al, J. Physiol. 441: 175, 1991.
:
: Kinetic model of calcium-induced shift in the activation of Ih channels
: Model of A. Destexhe, 1992, inspired from the dependence of If on calcium
: in heart cells (Harigawa & Hirishawa, J. Physiol. 409: 121, 1989)
:
:   ACTIVATE BINDING MODEL : 
:       - binding of Ca on S and F channels (VERSION 2: nexp binding sites)
:       - Ca binds on activated gates (rate constants k1 and k2)
:	    idem before:
:		s0 (closed) <-> s1 (open)	; rate cst alpha1,beta1
:		f0 (closed) <-> f1 (open)	; rate cst alpha1,beta1

:	    new:
:		s1 (open) + Ca <-> s2 (open)	; rate cst k1,k2
:		f1 (open) + Ca <-> f2 (open)	; rate cst k1,k2
:
:       - this suffies to account for shift of Ih activation with calcium
:	  (no need of other mechanism - or other time constants than k1,k2)
:
:   PARAMETERS:
:
:     VERSION 2: reformulation of parameters k1,k2 into k2 and cac.
:	cac = (k2/k1)^(1/nexp) = half activation calcium dependence.
:	- k2:  this rate constant is the inverse of the real time constant of 
:              the binding of Ca to Ih channel.  (0.001 to 0.0001 ms-1)
:	- cac: the half activation must be adapted to calcium dynamics of
:	       the cell.  Usually, cac = 1e-4 mM.
:	- nexp:number of sites of calcium on h-channels, nexp=2 here.
:
:  MODIF: addition of control variables (June 11 93)
:
: Written by Alain Destexhe, Salk Institute, Aug 1992
:
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX iar
	USEION other WRITE iother VALENCE 1
	USEION ca READ cai
        RANGE ghbar, gh, i
	GLOBAL k2, cac, nexp, h_inf, tau_s, tau_f, controls, controlf
}

UNITS {
	(molar)	= (1/liter)
	(mM)	= (millimolar)
	(mA) 	= (milliamp)
	(mV) 	= (millivolt)
	(msM)	= (ms mM)
}


PARAMETER {
	eh	= -43	(mV)
	celsius = 36	(degC)
	ghbar	= .0001	(mho/cm2)
	cac	= 1e-4	(mM)		: half-activation of calcium dependence
	k2	= 0.001	(1/ms)		: inverse of time constant
	nexp	= 2			: number of binding sites
	controls = 1			: control of variable s (0=no s1, s2)
	controlf = 1			: control of variable f (0=no f1, f2)
}


STATE {
	s1
	s2
	f1
	f2
}


ASSIGNED {
	v	(mV)
	cai	(mM)
	i	(mA/cm2)
	iother	(mA/cm2)
        gh	(mho/cm2)
	h_inf
	tau_s	(ms)
	tau_f	(ms)
	alpha1	(1/ms)
	alpha2	(1/ms)
	beta1	(1/ms)
	beta2	(1/ms)
	kk	(1/ms)
	fderiv	(1/ms)
	tadj
}


BREAKPOINT {
	SOLVE states METHOD runge

	if(controls == 0) {
		gh = ghbar * (f1+f2)
	} else if(controlf == 0) {
		gh = ghbar * (s1+s2)
	} else {
		gh = ghbar * (s1+s2) * (f1+f2)
	}
	
	i = gh * (v - eh)
	iother = i
}

DERIVATIVE states { LOCAL s0,f0
	evaluate_fct(v)

	s0 = 1 - s1 - s2
	f0 = 1 - f1 - f2

	kk = k2 * (5e-5/cac)^nexp

	fderiv = kk*s1 - k2*s2

	s1' = alpha1*s0 - beta1*s1 - fderiv
	s2' = fderiv

	fderiv = kk*f1 - k2*f2

	f1' = alpha2*f0 - beta2*f1 - fderiv
	f2' = fderiv
}

UNITSOFF
INITIAL {
:
:  Experiments of Coulter et al were at 36 deg.C
:  Q10 is assumed equal to 3
:
        tadj = 3.0 ^ ((celsius-36)/10)
	evaluate_fct(v)
	kk = k2 * (cai/cac)^nexp
	s1 = alpha1*k2/(alpha1*kk + alpha1*k2 + beta1*k2)
	s2 = alpha1*kk/(alpha1*kk + alpha1*k2 + beta1*k2)
	f1 = alpha2*k2/(alpha2*kk + alpha2*k2 + beta2*k2)
	f2 = alpha2*kk/(alpha2*kk + alpha2*k2 + beta2*k2)
}


PROCEDURE evaluate_fct(v (mV)) {

	h_inf = 1 / ( 1 + exp((v+68.9)/6.5) )	: sigmoide "square root"
	tau_s = exp((v+183.6)/15.24) / tadj	: version J neuro
	tau_f = exp((v+158.6)/11.2) / ( 1 + exp((v+75)/5.5) ) / tadj

	alpha1 = controls * h_inf / tau_s
	beta1  = ( 1 - h_inf ) / tau_s
	alpha2 = controlf * h_inf / tau_f
	beta2  = ( 1 - h_inf ) / tau_f
}
UNITSON

