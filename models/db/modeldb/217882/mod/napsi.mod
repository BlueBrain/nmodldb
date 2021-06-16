TITLE napsi		:modified to have slow inactivation described in Fleidervish and to make slope a global parameter
      			:further modified to have intermediate time course inactivation

NEURON {
	SUFFIX napsi
	USEION na READ ena WRITE ina
	RANGE  gbar, thegna, sh, ar,ari,minf,sinf,hinf,taui,taus
	GLOBAL mtau,vslope
}

PARAMETER {
	gbar = .0052085   	(mho/cm2)
	sh = 0  (mV)
	vslope=6.8   (mV)	:activation slope
	mtau = 1 (ms)
	ena		(mV)       :must be explicitly defined in hoc     
        a0s=0.001	(/ms)	
        b0s=0.0034	(/ms)
        asvh=-85	(mV) 
        bsvh=-17	(mV) 
        avs=30		(mV)
        bvs=10		(mV)
        ar=1		(1)		: 1=no inact., 0=max inact.
        a0i=0.01	(/ms)	
        b0i=0.034	(/ms)
        aivh=-85	(mV) 
        bivh=-17	(mV) 
        avi=30		(mV)
        bvi=10		(mV)
        ari=1		(1)		: 1=no inact., 0=max inact.
	celsius (degC)
	v 		(mV)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ina 		(mA/cm2)
	thegna		(mho/cm2)
	minf 	
	sinf
	hinf
	taus (ms)
	taui (ms)
}
 

STATE { m s h}

UNITSOFF

BREAKPOINT {
    SOLVE states METHOD cnexp
	        		
	thegna =gbar*m*s*h   
	ina = thegna * (v - ena)
	} 

INITIAL {
	trates(v,ar,ari,sh)
	m=minf  
	s=sinf
	h=hinf
}

DERIVATIVE states {   
    	trates(v,ar,ari,sh)
        s' = (sinf - s)/taus
	m' = (minf-m)/mtau
	h' = (hinf-h)/taui
}

PROCEDURE trates(vm,a1,a2,sh2) {  
        LOCAL   c,ci  

	minf = (1/(1+exp(-(vm+52.3-sh2)/vslope)))      	
        taus = 1/(alps(vm)+bets(vm))
        taui = 1/(alpi(vm)+beti(vm))
	c=alps(vm)*taus
        sinf = c+a1*(1-c)
	ci=alpi(vm)*taui
        hinf = ci+a2*(1-ci)
 }




FUNCTION alps(v(mV)) {  
  alps = a0s*exp((asvh-v)/avs)
}

FUNCTION bets(v(mV)) {
  bets = b0s/(exp((bsvh-v)/bvs)+1)
}

FUNCTION alpi(v(mV)) {  
  alpi = a0i*exp((aivh-v)/avi)
}

FUNCTION beti(v(mV)) {
  beti = b0i/(exp((bivh-v)/bvi)+1)
}

UNITSON
