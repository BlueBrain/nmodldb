TITLE I-h channel for Thalamic neurons from McCormick and Pape (1990)

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
        v 		(mV)
        eh  		(mV)        
        celsius 	(degC)
        ghbar=.00005 	(mho/cm2)
        vhalft=-80   	(mV)
        a0t=0.0005      	(/ms)
        zetat=0.2    	(1)
        gmt=.65   	(1)
        q10=4.5
        shift=0   (mV)
}


NEURON {
        SUFFIX htc
        NONSPECIFIC_CURRENT i
        RANGE ghbar, eh
	THREADSAFE
        GLOBAL linf,taul,shift
}

STATE {
        l
}

ASSIGNED {
        i (mA/cm2)
        linf      
        taul
}

INITIAL {
	rate(v)
	l=linf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	i = ghbar*l*(v-eh)

}


FUNCTION alpt(v(mV)) {
  alpt = exp(zetat*(v-vhalft+shift)) 
}

FUNCTION bett(v(mV)) {
  bett = exp(zetat*gmt*(v-vhalft+shift)) 
}

DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rate(v)
        l' =  (linf - l)/taul
}

PROCEDURE rate(v (mV)) { :callable from hoc
        LOCAL qt
        qt=q10^((celsius-36)/10)
        linf = 1/(1+ exp((v+75+shift)/5.5))
        taul = bett(v)/(qt*a0t*(1+alpt(v)))
}