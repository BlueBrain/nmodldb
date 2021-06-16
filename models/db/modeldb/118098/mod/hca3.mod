TITLE I-h channel from Magee 1998 for distal dendrites
: with tapering and 2-exp to model currents in CA3, M.Migliore Dec.2005

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v 		(mV)
        ehd  		(mV)        
	celsius 	(degC)
	ghdbar=.0001 	(mho/cm2)
        vhalfl=-82   	(mV)
	kl=-7.8
        vhalft=-65   	(mV)
        a0t=0.012      	(/ms)
        zetat=8    	(1)
        gmt=.15   	(1)
	q10=4.5
	qtl=1
	b0=20
	vc=-70
	kc=-3
	as=0.79
        vhalfts=-65   	(mV)
        a0ts=0.0019     	(/ms)
        zetats=8    	(1)
        gmts=.18   	(1)
	b0s=140
}


NEURON {
	SUFFIX hd
	NONSPECIFIC_CURRENT i
        RANGE ghdbar
        GLOBAL linf,taul, tauls
}

STATE {
        l
	ls
}

ASSIGNED {
	i (mA/cm2)
        linf      
        taul
        tauls
        ghd
}

INITIAL {
	rate(v)
	l=linf
	ls=linf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ghd = ghdbar*(ls*as+l*(1-as))
	i = ghd*(v-ehd)

}


FUNCTION alpt(v(mV)) {
  alpt = exp(0.0378*zetat*(v-vhalft)) 
}

FUNCTION bett(v(mV)) {
  bett = exp(0.0378*zetat*gmt*(v-vhalft)) 
}

FUNCTION alpts(v(mV)) {
  alpts = exp(0.0378*zetats*(v-vhalfts)) 
}

FUNCTION betts(v(mV)) {
  betts = exp(0.0378*zetats*gmts*(v-vhalfts)) 
}

DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rate(v)
        l' =  (linf - l)/taul
        ls' =  (linf - ls)/tauls
}

PROCEDURE rate(v (mV)) { :callable from hoc
        LOCAL a,qt
:        qt=q10^((celsius-33)/10)
        linf = (1/(1 + exp(-(v-vhalfl)/kl)))
        a = alpt(v)
        taul = b0 + bett(v)/(a0t*(1+a))
        a = alpts(v)
        tauls = b0s + betts(v)/(a0ts*(1+a))
}














