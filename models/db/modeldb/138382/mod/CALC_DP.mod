TITLE Cerebellum Golgi Cell Model

COMMENT
        Calcium first order kinetics
   
	Author: A. Fontana
	Revised: 12.12.98

Adapted by: Haroon Anwar (anwar@oist.jp)

ENDCOMMENT

NEURON {
        SUFFIX CALC_DP
        USEION ca READ ica, cao WRITE cai
	USEION ca2 READ ica2 VALENCE 2
        RANGE d, beta, cai0
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
        cao = 2.        (mM)         
        cai0 = 45e-6	(mM)         
        beta = 1.3      (/ms)
}

STATE {
	cai (mM)
}

INITIAL {
        cai = cai0 
}

BREAKPOINT {
       SOLVE conc METHOD derivimplicit
}

DERIVATIVE conc {    
	cai' =  -(ica+ica2)/(2*F*d)*(1e4) - beta*(cai-cai0)
}


