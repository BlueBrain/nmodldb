TITLE KIR channel

NEURON {
	SUFFIX KIR
	USEION  k READ ek WRITE ik
	RANGE g, ik
	GLOBAL minf, mtau, gmax
}

CONSTANT {
	Q10 = 3 (1)
}

UNITS {
	(mA) = (milliamp)
	(uA) = (microamp)
	(mV) = (millivolt)
	(mS) = (millimho)
}

PARAMETER {
	ek			(mV)
	gmax = 1.4e-4	(mho/cm2)	<0,1e9>
	m_vh = -82	(mV)	: half activation
	m_ve = 13		(mV)	: slope
}

ASSIGNED {
	v	(mV)
	g	(mho/cm2)
	ik	(mA/cm2)
	minf	(1)
	mtau	(ms)
	qt (1)
}

STATE {
	m
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gmax*m
	ik = g*(v - ek)
}

INITIAL {
	qt = Q10^((celsius-35)/10)
	rates(v)
	m = minf
}

DERIVATIVE states { 
	rates(v)
	m' = (minf-m)/mtau
}

FUNCTION_TABLE tabmtau(v(mV)) (ms)

: rates() computes rate and other constants at present v
: call once from hoc to initialize inf at resting v

PROCEDURE rates(v(mV)) {
:	mtau = tabmtau(v)
	mtau = tabmtau(v)/qt
	minf = 1/(1 + exp((v - m_vh)/m_ve))
}

