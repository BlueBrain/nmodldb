TITLE Cardiac IKs  current
: Hodgkin - Huxley type K channel, from Courtemanche et al Am J Physiol 1998 275:H301


NEURON {
	SUFFIX IKs
	USEION k READ ek WRITE ik
	RANGE gKs, ik
	GLOBAL minf, mtau 
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (mM) = (milli/liter)
}

PARAMETER {
	 gKs=0.258e-3 (S/cm2) <0,1e9>
	
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
	ik = gKs*m*m*(v - ek)
}

DERIVATIVE states {	: 
	rate(v*1(/mV))
	m' = (minf - m)/mtau
}

UNITSOFF
FUNCTION alp(v(mV)) { LOCAL q10 
	v = v
	q10 = 3^((celsius - 37)/10)
        alp = q10*4e-5*(v - 19.9)/(1 - exp(-(v - 19.9)/17))
          
}

FUNCTION bet(v(mV)) { LOCAL q10  
	v = v 
	q10 = 3^((celsius - 37)/10)
        bet = q10*3.5e-5*(v - 19.9)/( exp((v - 19.9)/9) - 1)
        
}
                
FUNCTION ce(v(mV)) { 
        v = v
      
       
         ce = 1/(1 + exp(-(v - 19.9)/12.7))^0.5
        
}


PROCEDURE rate(v) {LOCAL a,b,c :
	TABLE minf, mtau DEPEND celsius FROM -100 TO 100 WITH 200
		a = alp(v)  b = bet(v)  c = ce(v)
		mtau = 0.5/(a + b)
		minf = c
               
}
UNITSON
