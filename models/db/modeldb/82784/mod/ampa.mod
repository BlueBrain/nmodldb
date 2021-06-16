: from Durstewitz & Gabriel (2006), Cerebral Cortex

TITLE ampa synapse 

NEURON {
	POINT_PROCESS ampa
	NONSPECIFIC_CURRENT i
        RANGE g,a,b,gAMPAmax,tauD,tauF,util
}

UNITS {
        (uS) = (microsiemens)
        (nA) = (nanoamp)
        (mV) = (millivolt)
}

PARAMETER {
	tcon = .2 (ms)
	tcoff = 1.0 (ms)
	eampa = 0 	(mV)
	gAMPAmax = 0	(uS)
        tauD = 800         (ms)
        tauF = 800         (ms)
        util= .3
}

ASSIGNED {
	v 	(mV)
	i	(nA)
	g       (uS)
	factor
}

INITIAL { 
   a=0  
   b=0 
   factor=tcon*tcoff/(tcoff-tcon)
}

STATE {
      a
      b
}

BREAKPOINT {
	SOLVE states METHOD derivimplicit
        g = b-a
	i = gAMPAmax*g*(v-eampa)
}

DERIVATIVE states {
	a' = -a/tcon
	b' = -b/tcoff
}

NET_RECEIVE(wgt,R,u,tlast (ms),nspike) {
        LOCAL x
        if (nspike==0) { R=1  u=util }
	else {
	     if (tauF>0) { u=util+(1-util)*u*exp(-(t-tlast)/tauF) }
	     if (tauD>0) { R=1+(R*(1-u)-1)*exp(-(t-tlast)/tauD) }
	     }
	x=wgt*factor*R*u
	state_discontinuity(a,a+x)
	state_discontinuity(b,b+x)
        tlast=t
        nspike= nspike+1
}
