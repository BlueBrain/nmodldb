: kst_wustenberg.mod --- 
: 
: Filename: kst_wustenberg.mod
: Description: 
: Author: Subhasis Ray
: Maintainer: 
: Created: Wed Dec 13 19:06:03 EST 2017
: Version: 
: Last-Updated: Mon Jun 18 14:35:01 2018 (-0400)
:           By: Subhasis Ray
: URL: 
: Doc URL: 
: Keywords: 
: Compatibility: 
: 
: 

: Commentary: 
: 
: NEURON implementation of slow transient K+ channel ( KST ) from Wustenberg
:DG, Boytcheva M, Grunewald B, Byrne JH, Menzel R, Baxter DA

:This is slow transient K+ channel in Apis mellifera Kenyon cells
:(cultured).

TITLE Slow transient KST current in honey bee from Wustenberg et al 2004

COMMENT
  NEURON implementation by Subhasis Ray (ray dot subhasis at gmail dot com).
  This channel has the same activation and inactivation parameters
  as KA, but slower kinetics:

ENDCOMMENT

INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

NEURON { 
        SUFFIX kst
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
:IK,ST  -81     8.11    minf                                    -20.1   16.1                    3
:                       taum    5.0             0.5             20      20
:                       hinf                                    -74.7   7                       1
:                       tauh    200             150             52      15
:
: The equations are:
:       minf = 1 / ( 1 + exp((Vh - V) / s))
:       hinf = 1 / ( 1 + exp((V - Vh) / s))
:       taum = (taumax - taumin) / (1 + exp((V - Vh1) / s1)) * (1 + exp((V - Vh2) / s2)) + taumin

PROCEDURE settables(v (mV)) { 
UNITSOFF
        TABLE minf, hinf, mtau, htau FROM -120 TO 40 WITH 641
        minf  = 1.0 / (1 + exp((-20.1 - v)/16.1))
        hinf  = 1.0 / (1 + exp((v + 74.7) / 7 ))
        mtau = (5.0 - 0.5) / (1 + exp((v - 20) / 20.0)) + 0.5
        htau = (200.0 - 150.0) / (1 + exp((v - 52) / 15.0)) + 150
UNITSON
}





: 
: kst_wustenberg.mod ends here
