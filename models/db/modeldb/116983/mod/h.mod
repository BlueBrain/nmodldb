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
        vhalfl=-100   	(mV)
	kl=-12
        vhalft=-75   	(mV)
        a0ta=0.001      	(/ms)
        a0td=0.009      	(/ms)
        zetat=2.2    	(1)
        gmt=.4   	(1)
	q10=4.5
	qtl=1
}


NEURON {
	SUFFIX hd
	NONSPECIFIC_CURRENT i
        RANGE ghdbar, vhalfl, tau
        GLOBAL linf,taua, taud
}

STATE {
        l
}

ASSIGNED {
	i (mA/cm2)
        linf      
	tau
        taua
        taud
        ghd
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
        if (l<linf) {tau=taua} else {tau=taud}
        l' =  (linf - l)/tau
}

PROCEDURE rate(v (mV)) { :callable from hoc
        LOCAL a,qt
        qt=q10^((celsius-33)/10)
        a = alpt(v)
        linf = 1/(1 + exp(-(v-vhalfl)/kl))
        taua = bett(v)/(qtl*qt*a0ta*(1+a))
        taud = bett(v)/(qtl*qt*a0td*(1+a))
}














