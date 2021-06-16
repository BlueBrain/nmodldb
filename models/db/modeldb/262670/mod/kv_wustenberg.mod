: kv_wustenberg.mod --- 
: 
: Filename: kv_wustenberg.mod
: Description: 
: Author: Subhasis Ray
: Maintainer: 
: Created: Wed Dec 13 19:06:03 EST 2017
: Version: 
: Last-Updated: Mon Jun 18 14:36:04 2018 (-0400)
:           By: Subhasis Ray
: URL: 
: Doc URL: 
: Keywords: 
: Compatibility: 
: 
: 

: Commentary: 
: 
: NEURON implementation of A type K+ channel ( KV ) from Wustenberg
:DG, Boytcheva M, Grunewald B, Byrne JH, Menzel R, Baxter DA
:
: This isdelayed rectifier type K+ channel in Apis mellifera Kenyon cells
:(cultured).

TITLE Transient KV (delayed rectifier) current in honey bee from Wustenberg et al 2004

COMMENT
  NEURON implementation by Subhasis Ray (ray dot subhasis at gmail dot com)

ENDCOMMENT

INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

NEURON { 
        SUFFIX kv
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
        ek      (mV)
        v       (mV)
        ik      (mA/cm2)
        g       (S/cm2)
        minf
        mtau    (ms)
}
 
STATE {
    m
}

BREAKPOINT { 
        SOLVE states METHOD cnexp 
        g = gbar * m * m * m * m
        ik = g * ( v - ek ) 
}
 
INITIAL { 
        settables(v)
        m = minf
} 

DERIVATIVE states { 
        settables(v) 
        m' = (minf - m ) / mtau
}

: Parameters from the article (Table 2):
:
:       E,mV    g,nS            taumax,ms       taumin,ms       Vh1,mV  s1      Vh2,mV  s2      N
:IK,V   -81     6       minf                                    -37.6   27.24                   4
:                       taum    3.53            1.85             45     13.71
: The equations are:
:       minf = 1 / ( 1 + exp((Vh - V) / s))
:       hinf = 1 / ( 1 + exp((V - Vh) / s))
:       taum = (taumax - taumin) / (1 + exp((V - Vh1) / s1)) + taumin

PROCEDURE settables(v (mV)) { 
UNITSOFF
        TABLE minf, mtau FROM -120 TO 40 WITH 641
        minf  = 1.0 / (1 + exp((-37.6 - v) / 27.24))
        mtau = (3.53 - 1.85) / (1 + exp((v - 45) / 13.71)) + 1.85
UNITSON
}





: 
: kv_wustenberg.mod ends here
