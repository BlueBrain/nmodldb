TITLE BK Ca 2+ -activated K + channel
: Calcium activated K channel.
COMMENT
 Starting from the formulation in De Schutter and Bower, 1994, we
reduced the Ca 2+ dependent activation time to half to account for the larger slow repolarisation at
depolarised states.
Current Model Reference: Karima Ait Ouares , Luiza Filipis , Alexandra Tzilivaki , Panayiota Poirazi , Marco Canepari (2018) Two distinct sets of Ca 2+ and K + channels 
are activated at different membrane potential by the climbing fibre synaptic potential in Purkinje neuron dendrites. 

PubMed link: 

Contact: Filipis Luiza (luiza.filipis@univ-grenoble-alpes.fr)

ENDCOMMENT

UNITS {
	(molar) = (1/liter)
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX KC3
	USEION ca READ cai
	USEION k WRITE ik
	RANGE gkbar,gk,zinf,ik,tin
}


PARAMETER {
	celsius=37	(degC)
	v		(mV)
	gkbar=.08	(mho/cm2)	: Maximum Permeability
	cai = .04e-3	(mM)
	ek  = -85	(mV)
	dt		(ms)
	tin=4
}


ASSIGNED {
	ik		(mA/cm2)
	minf
	mexp
	zinf
	zexp
	gk
}

STATE {	m z }		: fraction of open channels

BREAKPOINT {
	SOLVE state
:	gk = gkbar*1000*m*z*z
	ik = gkbar*1000*m*z*z*(v - ek)
}
:UNITSOFF
:LOCAL fac

:if state_cagk is called from hoc, garbage or segmentation violation will
:result because range variables won't have correct pointer.  This is because
: only BREAKPOINT sets up the correct pointers to range variables.
PROCEDURE state() {	: exact when v held constant; integrates over dt step
	rate(v, cai)
	m = m + mexp*(minf - m)
	z = z + zexp*(zinf - z)
	VERBATIM
	return 0;
	ENDVERBATIM
}

INITIAL {
	rate(v, cai)
	m = minf
	z = zinf
}

FUNCTION alp(v (mV), ca (mM)) (1/ms) { :callable from hoc
	alp = 400/(ca*1000)
}

FUNCTION bet(v (mV)) (1/ms) { :callable from hoc
	bet = 0.11/exp((v-55)/14.9)
}

PROCEDURE rate(v (mV), ca (mM)) { :callable from hoc
	LOCAL a,b
	a = alp(v,ca)
	zinf = 1/(1+a)
	zexp = (1 - exp(-dt/tin))
	b = bet(v)
	minf = 8.5/(7.5+b)
	mexp = (1 - exp(-dt*(7.5+b)))
}
:UNITSON
