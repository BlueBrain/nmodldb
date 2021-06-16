TITLE CA1 KM channel from M. Taglialatela, Kv72wt+Kv72R201C+Kv73wt
: M. Migliore Oct 2014

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v 		(mV)
	ek
	celsius 	(degC)
	gbar=.0001 	(mho/cm2)
        vhalfl=-43.4   	(mV)
		kl=-13.3
        vhalft=-67   	(mV)
        a0a=0.0035      	(/ms)
        zetat=5   	(1)
        gmt=.66   	(1)
        vhalfb=-73   	(mV)
        a0b=0.012      	(/ms)
        zetab=5    	(1)
        gmb=.9   	(1)
	q10=3.8
	b0=106
	b0b=30
	}


NEURON {
	SUFFIX kv72wt73wtR201C
	USEION k READ ek WRITE ik
        RANGE  gbar,ik
      GLOBAL inf, tau, taua, taub
}

STATE {
        m
}

ASSIGNED {
	ik (mA/cm2)
        inf
	tau
    taua
	taub
}

INITIAL {
	rate(v)
	m=inf
}


BREAKPOINT {
	SOLVE state METHOD cnexp
	ik = gbar*m*(v-ek)
}


FUNCTION alpa(v(mV)) {
  alpa = exp(0.0378*zetat*(v-vhalft)) 
}

FUNCTION alpb(v(mV)) {
  alpb = exp(0.0378*zetab*(v-vhalfb)) 
}


FUNCTION beta(v(mV)) {
  beta = exp(0.0378*zetat*gmt*(v-vhalft)) 
}

FUNCTION betb(v(mV)) {
  betb = exp(0.0378*zetab*gmb*(v-vhalfb)) 
}


DERIVATIVE state {
    rate(v)
    if (m<inf) {tau=taua} else {tau=taub}
	m' = (inf - m)/tau
}

PROCEDURE rate(v (mV)) { :callable from hoc
        LOCAL a,qt, ab, ac
        qt=q10^((celsius-22)/10)
        inf = (1/(1 + exp((v-vhalfl)/kl)))
        a = alpa(v)
        ab = alpb(v)
        taua = (b0 + beta(v)/(a0a*(1+a)))/qt
        taub = (b0b + betb(v)/(a0b*(1+ab)))/qt
}














