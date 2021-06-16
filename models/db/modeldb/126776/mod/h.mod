TITLE I-h channel from Magee 1998 for distal dendrites
: modified to take into account Sonia's exp. Apr.2008 M.Migliore

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
        vhalft=-62   	(mV)
        a0t=0.0077696      	(/ms)
        zetat=5    	(1)
        gmt=.057127   	(1)
	q10=4.5
	qtl=1
}


NEURON {
	SUFFIX hd
	NONSPECIFIC_CURRENT i
        RANGE ghdbar, vhalfl
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
        qt=q10^((celsius-33)/10)
        a = alpt(v)
        linf = 1/(1 + exp(-(v-vhalfl)/kl))
:       linf = 1/(1+ alpl(v))
        taul = bett(v)/(qtl*qt*a0t*(1+a))
}














