TITLE Cardiac IKr current
: Hodgkin - Huxley type K channel, from Courtemanche et al Am J Physiol 1998 275:H301


NEURON {
	SUFFIX IKr
	USEION k READ ek WRITE ik
	RANGE gKr, ik, Tauact
	GLOBAL minf, mtau 
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (mM) = (milli/liter)
	
}

PARAMETER {
	 gKr=0.0588e-3 (S/cm2) <0,1e9>
	Tauact=1 (ms)
	
}

STATE {
	 m 
}

ASSIGNED {
	v (mV)
	celsius (degC) : 37
	ik (mA/cm2)
	minf 
	mtau (ms)  
	ek (mV)      
}

INITIAL {
	rate(v*1(/mV))
	m = minf
}

BREAKPOINT {
	SOLVE states METHOD derivimplicit
	ik = gKr/(1 + exp((v + 15)/22.4))*m*(v - ek)
}

DERIVATIVE states {	
	rate(v*1(/mV))
	m' = (minf - m)/mtau
}

UNITSOFF
FUNCTION alp(v(mV)) { LOCAL q10 
	v = v
	q10 = 3^((celsius - 37)/10)
        alp = q10*0.0003*(v + 14.1)/(1 - exp(-(v + 14.1)/5))
          
}

FUNCTION bet(v(mV)) { LOCAL q10  
	v = v 
	q10 = 3^((celsius - 37)/10)
        bet = q10*7.3898e-5*(v - 3.3328)/( exp((v - 3.3328)/5.1237) - 1)
        
}
                
FUNCTION ce(v(mV)) {  
        v = v
       
       
         ce = 1/(1 + exp(-(v + 14.1)/6.5))
        
}


PROCEDURE rate(v) {LOCAL a,b,c :
	:TABLE minf, mtau DEPEND celsius FROM -100 TO 100 WITH 200
		a = alp(v)  b = bet(v)  c = ce(v)
		mtau = 1/(a + b)*Tauact
		minf = c
               
}
UNITSON
