TITLE I-h channel from Magee 1998 for distal dendrites

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v 		(mV)
    ehd  		(mV)        
	celsius 	(degC)
	ghdbar=.0001 	(mho/cm2)
    vhalfl=-81   	(mV)
	kl=-8
    vhalft=-75   	(mV)
    a0t=0.011      	(/ms)
    zetat=2.2    	(1)
    gmt=.4   	(1)
	q10=4.5
	qtl=1
	tfactor=1	: For running simulations with various peak values of tau_h
	consttau=50	(ms) : For running simulations with constant tau_h
}


NEURON {
	SUFFIX hd
	NONSPECIFIC_CURRENT i
        RANGE ghdbar, vhalfl,tfactor, ghd, ih, consttau
        GLOBAL linf,taul
}

STATE {
        l
}

ASSIGNED {
	i (mA/cm2)
        linf      
        taul
        ghd
		ih
}

INITIAL {
	rate(v)
	l=linf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ghd = ghdbar*l
	i = ghd*(v-ehd)
	ih = ghd*(v-ehd)
}


FUNCTION alpt(v(mV)) {
  alpt = exp(0.0378*zetat*(v-vhalft)) 
}

FUNCTION bett(v(mV)) {
  bett = exp(0.0378*zetat*gmt*(v-vhalft)) 
}

DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rate(v)
        l' =  (linf - l)/taul
}

PROCEDURE rate(v (mV)) { :callable from hoc
        LOCAL a,qt
        qt=q10^((celsius-33)/10)
        a = alpt(v)
        linf = 1/(1 + exp(-(v-vhalfl)/kl))
:       linf = 1/(1+ alpl(v))
 		taul = tfactor*bett(v)/(qtl*qt*a0t*(1+a))

: If you want to use consttau, comment the line above and uncomment the
: one below, and assign/change consttau in the HOC files.

:		taul = consttau
}
