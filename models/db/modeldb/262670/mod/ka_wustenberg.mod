: ka_wustenberg.mod --- 
: 
: Filename: ka_wustenberg.mod
: Description: 
: Author: Subhasis Ray
: Maintainer: 
: Created: Wed Dec 13 19:06:03 EST 2017
: Version: 
: Last-Updated: Mon Jun 18 14:33:52 2018 (-0400)
:           By: Subhasis Ray
: URL: 
: Doc URL: 
: Keywords: 
: Compatibility: 
: 
: 

: Commentary: 
: 
: NEURON implementation of A type K+ channel ( KA ) from Wustenberg
: DG, Boytcheva M, Grunewald B, Byrne JH, Menzel R, Baxter DA This is
: transient A type K+ channel in Apis mellifera Kenyon cells
: (cultured).

TITLE Transient KA current in honey bee from Wustenberg et al 2004

COMMENT
  NEURON implementation by Subhasis Ray (ray dot subhasis at gmail dot com)

ENDCOMMENT

INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

NEURON { 
        SUFFIX ka
        USEION k READ ek WRITE ik
        RANGE gbar, ik, g
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
	ek	(mV)
        v	(mV)
        ik	(mA/cm2)
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
        ik = g * ( v - ek )
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
:IK,A   -81     58.1    minf                                    -20.1   16.1                    3
:                       taum    1.65             0.35           -70     -4      -20     12      
:                       hinf                                    -74.7   7                       1
:                       tauh    90               2.5            -60     -25     -62     16
:
: The equations are:
:       minf = 1 / ( 1 + exp((Vh - V) / s))
:       hinf = 1 / ( 1 + exp((V - Vh) / s))
:       taum = (taumax - taumin) / (1 + exp((V - Vh1) / s1)) * (1 + exp((V - Vh2) / s2)) + taumin

PROCEDURE settables(v (mV)) { 
UNITSOFF
        TABLE minf, hinf, mtau, htau FROM -120 TO 40 WITH 641
        minf  = 1.0 / (1 + exp((-20.1 - v)/16.1))
        hinf  = 1.0 / ( 1 + exp( ( v + 74.7 ) / 7 ) )
        mtau = (1.65 - 0.35) / ((1 + exp(- (v + 70) / 4.0)) * (1 + exp((v + 20) / 12.0))) + 0.35
        htau = (90 - 2.5) / ((1 + exp(- (v + 60) / 25.0)) * (1 + exp((v + 62) / 16.0))) + 2.5
UNITSON
}





: 
: ka_wustenberg.mod ends here
