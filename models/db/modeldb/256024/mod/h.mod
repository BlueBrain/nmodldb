TITLE I-h channel from Magee 1998 for distal dendrites

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	ghdbar=1e-4	(mho/cm2)
	vhalfl=-81	(mV)
	kl=-8		(mV)
	vhalft=-75	(mV)
	a0t=0.011	(/ms)
	zetat=2.2	(/mV)
	gmt=.4   	(1)
	q10=4.5		(1)
	qtl=1		(1)
	ehd=-30  	(mV)        

}


NEURON {
    THREADSAFE
    
	SUFFIX hd
	NONSPECIFIC_CURRENT i
	RANGE ghdbar, ghd, taul, vhalfl
	:GLOBAL linf,taul
}

STATE {
        l
}

ASSIGNED {
	v 		(mV)
	celsius (degC)
	i		(mA/cm2)
	linf	(1)
	taul	(ms)
	ghd		(mho/cm2)
}

INITIAL {
	rate(v)
	l=linf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ghd = ghdbar*l
	i = ghd*(v-ehd)

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
        qt=q10^((celsius-33(degC))/10(degC))
        a = alpt(v)
        linf = 1/(1 + exp(-(v-vhalfl)/kl))
:       linf = 1/(1+ alpl(v))
        taul = bett(v)/(qtl*qt*a0t*(1+a))
}














