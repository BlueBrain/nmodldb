TITLE CAPUMP
 

UNITS {
       (molar) = (1/liter)
        (pA) = (picoamp)
	(mV) =	(millivolt)
        (uS) = (micromho)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}


INDEPENDENT {v FROM -100 TO 50 WITH 50 (mV)}

NEURON {
	SUFFIX capump
	USEION ca READ cai  WRITE ica
	RANGE  icapump,icapumpmax,n,km
 
}


PARAMETER {
        dt (ms)
        cai   (mM)
        icapumpmax  = 0.0312   (mA/cm2)
        km = 0.000500    (mM)
        n=1.0
        cainit = 0.000184 (mM)
        celsius = 35  (degC)
        
 
}

ASSIGNED { 
           ica		(mA/cm2)
        icapump (mA/cm2)
}


BREAKPOINT {
        icapump = icapumpmax*(1/(1 + km/cai))
	ica = icapump
}


COMMENT
INITIAL{
       cai = cainit}
ENDCOMMENT
