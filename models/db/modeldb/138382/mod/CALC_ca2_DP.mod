TITLE Cerebellum Golgi Cell Model

COMMENT
        Calcium first order kinetics
   
	Author: A. Fontana
	Revised: 12.12.98

Adapted by: Haroon Anwar (anwar@oist.jp)

ENDCOMMENT

NEURON {
        SUFFIX CALC2
	USEION ca READ ica
        USEION ca2 READ ica2 WRITE ca2i VALENCE 2
        RANGE d, beta, ca2i0 
}

UNITS {
        (mV)    = (millivolt)
        (mA)    = (milliamp)
	(um)    = (micron)
	(molar) = (1/liter)
        (mM)    = (millimolar)
   	F      = (faraday) (coulomb)
}

PARAMETER {
        ica             (mA/cm2)
	ica2		(mA/cm2)
        celsius    	(degC)
        d = .2          (um)
        ca2i0     	(mM)         
        beta = 1.3      (/ms)
}

STATE {
	ca2i (mM)
}

INITIAL {
        ca2i = ca2i0 
}

BREAKPOINT {
       SOLVE conc METHOD derivimplicit
}

DERIVATIVE conc {    
	ca2i' = -(ica+ica2)/(2*F*d)*(1e4) - beta*(ca2i-ca2i0)
}


