TITLE Cardiac fast sodium current
COMMENT
	modified From DiFrancesco & Noble 1985 Phil Trans R Soc Lond 307:353-398 
    modified for Neuron by FE GANNIER
	francois.gannier@univ-tours.fr (University of TOURS)
ENDCOMMENT
INCLUDE "Unit.inc"
INCLUDE "Volume.inc"
NEURON {
	SUFFIX INaf
	USEION na READ nao, nai, ena WRITE ina
	USEION k READ ko, ki

	RANGE Emh, g, ina, m_a, m_b, h_a, h_b, minf, mtau, hinf, htau
}

PARAMETER {
	g = 750		(uS)
}

STATE { : m, h
	m h
}

ASSIGNED {
	v (mV)
	celsius (degC) : 37
	ina (mA/cm2)
	ena (mV)
	Emh (mV)
	minf hinf
	mtau (ms)
	htau (ms)
	
	m_a (/ms)
	m_b (/ms)
	h_a (/ms)
	h_b  (/ms)
	
	ko (mM)
	ki (mM)
	nai (mM)
	nao (mM)
}

LOCAL RT
INITIAL {
	RT = (1000)*R*(273.15+celsius)
	rate(v)
	m = minf
	h = hinf
}

BREAKPOINT { 
	SOLVE states METHOD derivimplicit
	Emh = (RT/F) * log((nao + 0.12 * ko)/(nai + 0.12 * ki))
	ina = (1e-06)* g/S *m*m*m*h *(v - Emh)
}

DERIVATIVE states {
	rate(v)
:	m' = (minf - m)/mtau
	m' = m_a*(1 - m) - m_b*m
:	h' = (hinf - h)/htau
	h' = h_a*(1 - h) - h_b * h
}

FUNCTION alp(v(mV), i) (/ms) { LOCAL E0 : order m,h
	if (i==0) {
		E0 = (v + 41)
		if (fabs(E0*1(/mV)) < 1e-5)
		{
			alp = (0.001)* 2000 (/s)
		}
		else
		{
			alp = (0.001)*200(/s/mV)*E0 / (1 - exp(-0.1(/mV) * E0))
		}
	} else if (i==1) {
		alp = (0.001)*20(/s)*exp(-0.125(/mV)*(v + 75))
	}
}

FUNCTION bet(v(mV), i) (/ms) { : order m,h
	if (i==0) {
		bet = (0.001)*8000(/s) * exp(-0.056(/mV)*(v + 66))
	} else if (i==1) {
:		bet = 2000/(320*exp(-0.1*(v + 75)))
: correction
		bet = (0.001)*2000(/s)/(320*exp(-0.1(/mV)*(v + 75))+1)
	}
}

PROCEDURE rate(v(mV)) { 
	: TABLE minf, mtau DEPEND celsius FROM -100 TO 100 WITH 200
	m_a = alp(v,0)  m_b = bet(v,0)
	mtau = 1/(m_a + m_b)
	minf = m_a * mtau
	h_a = alp(v,1)  h_b = bet(v,1)
	htau = 1/(h_a + h_b)
	hinf = h_a * htau
}
