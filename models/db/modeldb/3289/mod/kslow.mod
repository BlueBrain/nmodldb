TITLE K-slow channel from Korngreen and Sakmann 2000
: M.Migliore June 2001

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v (mV)
	celsius		(degC)
	gkbar=.008 (mho/cm2)
        vhalfn=-14   (mV)
        vhalfl=-54   (mV)
        kn=14.6   (1)
        kl=-11   (1)
	q10=2.3
	tmp=13
	ek
}


NEURON {
	SUFFIX ks
	USEION k READ ek WRITE ik
        RANGE gkbar
        GLOBAL ninf,linf,taul,taun
}

STATE {
	n
        l
}

ASSIGNED {
	ik (mA/cm2)
        ninf
        linf      
        taul
        taun
}

INITIAL {
	rates(v)
	n=ninf
	l=linf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkbar*n^2*l*(v-ek)
}


DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rates(v)
        n' = (ninf - n)/taun
        l' =  (linf - l)/taul
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a,qt
        qt=q10^((celsius-22)/10)
        ninf = 1/(1 + exp(-(v-vhalfn)/kn))
        linf = 1/(1 + exp(-(v-vhalfl)/kl))
	if (v<-50) {taun = (1.25+175.03*exp(v*0.026))/qt
		   } else {
		    taun = (1.25+13*exp(-v*0.026))/qt}
        taul = (360+(1010+24*(v+55))*exp(-((v+75)/48)^2))/qt
}














