TITLE C channel
: Ca-dependent K channel.
: By Ojvind Bernander 92-01-21.
: (db) added RANGE section to allow access to parameters from mech. browser
: (db) 4.1.98 modifications for CVode

NEURON {
	SUFFIX icnew
	USEION k READ ek WRITE ik
	USEION ca READ ica
	GLOBAL inf
	RANGE gkbar,ek,tau_diff,taum,fact, ik, ca_beta
}
UNITSON
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
PARAMETER {
	v (mV)
	dt (ms)
	gkbar=.045 (mho/cm2)
	ek = -95 (mV)
	:ca_beta  = 20.0 (1/ms)   : Ca decay (inverse, 50 msec)
        tau_diff = 0.05 (ms)
	taum = 2.0 (1/ms)
	fact = 1e-3

}

CONSTANT {
	ca_beta  = 20.0 (1/ms)   : Ca decay (inverse, 50 msec)
	ca_alpha = 100.0 (mM/ms/mA)
}

STATE { m cai}
ASSIGNED {
	ik (mA/cm2)
	ica (mA/cm2)
	inf[1]
	tau[1]
}

INITIAL {
	cai = 0.0 (mM)
  	:ca_beta=1/(tau_diff+1e-10)
	:printf("ca_beta=%g\n",ca_beta)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkbar*m*m*(v - ek)
}

DERIVATIVE states {	: exact when ca held constant
	mhn(cai)
	m' = (inf[0] - m)/tau[0]
	cai' = (ca_alpha*(-ica) - ca_beta*cai)*fact
}

FUNCTION varss(ca) {
	varss = ca / (ca + 0.040) :K activation
	:printf("ca_beta=%g\n",ca_beta)
}

FUNCTION vartau() {
	vartau = taum  :2.0  K activation tau
}	


PROCEDURE mhn(ca) { :no dependence on voltage, only on ca concentration (appears in states)
:	TABLE inf,tau 
:	DEPEND celsius, dt
:	FROM -100 TO 100 WITH 2000 : .1 mV steps
	tau[0] = vartau()
	inf[0] = varss(ca)
}
