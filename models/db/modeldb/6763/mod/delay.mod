TITLE LEAK
 

UNITS {
        (pA) = (picoamp)
        (molar) = (1/liter)
	(mV) =	(millivolt)
        (uS) = (micromho)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
        F = (faraday) (coulomb)
        R = (mole k) (mV-coulomb/degC)
}


INDEPENDENT {v FROM -100 TO 50 WITH 50 (mV)}

NEURON {
	SUFFIX delay
	USEION k  WRITE ik
	RANGE  ik,gkbar,slope,vhalf
        GLOBAL ninf
 
}


PARAMETER {
        dt (ms)
        gkbar = 800  (uS/cm2)
        ek =  -100 (mV)
        celsius = 35  (degC)
        vhalf = 10 (mV) 
        slope = 10 (mV)
        
 
}

ASSIGNED { 
           ik		(mA/cm2)
        ninf
}


BREAKPOINT {
        ninf= 1/(1+ exp((vhalf-v)/slope))
	ik = (0.000001)*ninf*gkbar*(v-ek)
}


COMMENT
INITIAL{
       nai = nainit}
        ena = R*(celsius+273.15)/F*log(nao/nai)
ENDCOMMENT
