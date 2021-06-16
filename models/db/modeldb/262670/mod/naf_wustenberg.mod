: naf_wustenberg.mod --- 
: 
: Filename: naf_wustenberg.mod
: Description: 
: Author: Subhasis Ray
: Maintainer: 
: Created: Wed Dec 13 19:06:03 EST 2017
: Version: 
: Last-Updated: Mon Jun 18 14:37:06 2018 (-0400)
:           By: Subhasis Ray
: URL: 
: Doc URL: 
: Keywords: 
: Compatibility: 
: 
: 

: Commentary: 
: 
: NEURON implementation of fast Na+ channel ( NAF ) from Wustenberg
: DG, Boytcheva M, Grunewald B, Byrne JH, Menzel R, Baxter DA

: This is fast Na+ channel in Apis mellifera Kenyon cells :(cultured).

TITLE Fast NA+ current in honey bee from Wustenberg et al 2004

COMMENT
  NEURON implementation by Subhasis Ray (ray dot subhasis at gmail dot com).

ENDCOMMENT

INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

NEURON { 
        SUFFIX naf
        USEION na READ ena WRITE ina
        RANGE gbar, ina, g
}

UNITS {
        (S) = (siemens)
        (mV) = (millivolt) 
        (mA) = (milliamp) 
}
 
PARAMETER { 
        gbar = 0.0      (mho/cm2)
}
 
ASSIGNED { 
	ena	(mV)
        v	(mV)
        ina	(mA/cm2)
        g	(S/cm2)
        minf
	hinf
        mtau	(ms)
        htau	(ms)
}
 
STATE {
    m
    h
}

BREAKPOINT { 
        SOLVE states METHOD cnexp 
        g = gbar * m * m * m * h
        ina = g * ( v - ena )
}
 
INITIAL { 
        settables(v)
	m = minf
        h  = hinf
} 

DERIVATIVE states { 
        settables(v) 
        h' = (hinf - h) / htau
	m' = (minf - m ) / mtau
}

: Parameters from the article (Table 2):
:
:       E,mV    g,nS            taumax,ms       taumin,ms       Vh1,mV  s1      Vh2,mV  s2      N
: INa                     minf                                  -30.1   6.65                    3
:                         taum  0.83            0.093           -20.3   6.45
: INaF  58      140       hinf                                  -51.4   5.9                     1
:                         tauh  1.66            0.12            -8.03   8.69
:
: The equations are:
:       minf = 1 / ( 1 + exp((Vh - V) / s))
:       hinf = 1 / ( 1 + exp((V - Vh) / s))
:       taum = (taumax - taumin) / (1 + exp((V - Vh1) / s1)) + taumin

PROCEDURE settables(v (mV)) { 
UNITSOFF
        TABLE minf, hinf, mtau, htau FROM -120 TO 40 WITH 641
        minf  = 1.0 / (1 + exp((-30.1 - v) / 6.65))
        hinf  = 1.0 / (1 + exp((v + 51.4) / 5.9 ))
        mtau = (0.83 - 0.093) / (1 + exp((v + 20.3) / 6.45)) + 0.093
        htau = (1.66 - 0.12) / (1 + exp((v + 8.03) / 8.69)) + 0.12
UNITSON
}





: 
: naf_wustenberg.mod ends here
