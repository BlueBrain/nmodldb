: nattxs.mod is a transient ttx-sensitive Na+ current from
: Herzog et al. 2001 p1353
: This current likely consists of NaV1.7 with a little
: NaV1.6 and NaV1.1 mixed in (TMM).
: implemented by Tom Morse version 2/25/07

NEURON {
	SUFFIX nattxs
	NONSPECIFIC_CURRENT i
	RANGE gbar, ena, m, h, tau_m, tau_h
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 0.035135 (S/cm2)
	ena=62.94 (mV)

: second commented values are those used in Baker '05
	A_am = 11.49 (/ms)  : 17.235 (/ms) : A for alpha m
	B_am = 8.58 (mV)    : 7.58 (mV)
	C_am = -8.47 (mV)   : -11.47 (mV)

	A_ah = 0.0658 (/ms) : 0.23688 (/ms) : A for alpha h
	B_ah = 120 (mV)     : 115 (mV)
	C_ah = 20.33 (mV)   : 46.33 (mV)

	A_bm = 11.49 (mV)   : 17.235 (/ms) : A for beta m
	B_bm = 67.2 (mV)    : 66.2 (mV)
	C_bm = 27.8 (mV)    : 19.8 (mV)

	A_bh = 3.0 (/ms)    : 10.8 (/ms)   : A for beta h
	B_bh = -6.8 (mV)    : -11.8 (mV)
	C_bh = -12.998 (mV) : -11.998 (mV)
}

ASSIGNED {
	v	(mV) : NEURON provides this
	i	(mA/cm2)
	g	(S/cm2)
	tau_h	(ms)
	tau_m	(ms)
	minf
	hinf
}

STATE { m h }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * m^3 * h
	i = g * (v-ena)
}

INITIAL {
	rates(v) : set tau_m, tau_h, hinf, minf
	: assume that equilibrium has been reached
	m = minf
	h = hinf
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/tau_m
	h' = (hinf - h)/tau_h
}

FUNCTION alpham(Vm (mV)) (/ms) {
	alpham=A_am/(1+exp((Vm+B_am)/C_am))
}

FUNCTION alphah(Vm (mV)) (/ms) {
	alphah=A_ah*exp(-(Vm+B_ah)/C_ah)
}

FUNCTION betam(Vm (mV)) (/ms) {
	betam=A_bm/(1+exp((Vm+B_bm)/C_bm))
}

FUNCTION betah(Vm (mV)) (/ms) {
	betah=A_bh/(1+exp((Vm+B_bh)/C_bh))
}

FUNCTION rates(Vm (mV)) (/ms) {
	tau_m = 1.0 / (alpham(Vm) + betam(Vm))
	minf = alpham(Vm) * tau_m

	tau_h = 1.0 / (alphah(Vm) + betah(Vm))
	hinf = alphah(Vm) * tau_h
}
