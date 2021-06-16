: kdr.mod is the delayed rectifier K+ current from
: Herzoz, Cummins, and Waxman '01, parameter assignments and formula's
: from page 1353
: implemented by Tom Morse version 2/25/07

NEURON {
	SUFFIX kdr
	NONSPECIFIC_CURRENT i
	RANGE gbar, ek
	RANGE tau_n, n
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 0.0021 (S/cm2)
	: ek=-70 (mV) : this does not represent value in paper
      ek = -92.34 (mV) : correction by Tom Andersson
: Baker 2005 values
	A_anF = 0.001265 (/ms) : 0.00798 Baker '05 : A for alpha n
	B_anF = 14.273 (mV) : 72.2 Baker '05
	C_anF = 10 (mV) : 1.1 Baker '05

: Baker '05 uses different (parameterized) function than below beta_n
	A_bnF = 0.125 (/ms): A for beta n
	B_bnF = 55 (mV)
	C_bnF = -2.5 (mV)

: Bostok et al. 1991 values
:	A_anF = 0.129 (/ms) : A for alpha n
:	B_anF = -53 (mV)
:	C_anF = 10 (mV)

:	A_bnF = 0.324 (/ms) : A for beta n
:	B_bnF = -78 (mV)
:	C_bnF = 10 (mV)
}

ASSIGNED {
	v	(mV) : NEURON provides this
	i	(mA/cm2)
	g	(S/cm2)
	tau_n	(ms)
	ninf
}

STATE { n }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * n
	i = g * (v-ek)
}

INITIAL {
	: assume that equilibrium has been reached
	n = alphan(v)/(alphan(v)+betan(v))
}

DERIVATIVE states {
	rates(v)
	n' = (ninf - n)/tau_n
}

FUNCTION alphan(Vm (mV)) (/ms) {
	if (-Vm-B_anF != 0) {
		alphan=A_anF*(Vm+B_anF)/(1-exp((-Vm-B_anF)/C_anF))
	} else {
		alphan=A_anF*C_anF
	}
}

FUNCTION betan(Vm (mV)) (/ms) {
	betan=A_bnF*exp((Vm+B_bnF)/C_bnF)
}

FUNCTION rates(Vm (mV)) (/ms) {
	tau_n = 1.0 / (alphan(Vm) + betan(Vm))
:	ninf = alphan(Vm) * tau_n : this line does not reflect p. 1353
      ninf = 1/(1+exp((Vm+14.62)/-18.38)) : correction by Tom Andersson
}
