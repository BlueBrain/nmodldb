TITLE K-fast channel from Korngreen and Sakmann 2000
: M.Migliore June 2001

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v (mV)
	celsius		(degC)
	gkbar=.008 (mho/cm2)
        vhalfn=-47   (mV)
        vhalfl=-66   (mV)
        kn=29   (1)
        kl=-10   (1)
	q10=2.3
	qq=5
	tq=-55
	ek
}


NEURON {
	SUFFIX kf
	USEION k READ ek WRITE ik
        RANGE gkbar
        GLOBAL ninf,linf,taul,taun, tq,qq, q10
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
	ik = gkbar*n^4*l*(v-ek)
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
        taun = (0.34+0.92*exp(-((v+71)/59)^2))/qt
        taul = (8+49*exp(-((v+73)/23)^2))/qt
}














