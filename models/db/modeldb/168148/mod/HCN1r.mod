TITLE HCN1 with one time constant 
    :Konstantin Stadler 2009    
	
	
	
	
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {	
	v		(mV)			
	Vrev  = -40	(mV)
	gpeak = 0.00015	(mho/cm2) <0, 1e9>
		
	vhakt = -67.124	(mV)
	k     = -9.399	(mV)
		
	vhtau = -57	(mV)
	a0t   = 0.00682		(/ms)
	zetat = 1.8		(/mV)
	gmt   = .44		(1)
	
	celsius		(degC)			
	temp  = 24	(degC)
	q10   = 4.5		(1)
	qtl   = 1		(1)
}			


NEURON {
	SUFFIX hcn1
	NONSPECIFIC_CURRENT i
    	RANGE gpeak, vhakt,ghd		
}					
					
STATE {
        l
}

ASSIGNED {
	i (mA/cm2)
    	linf (1)     
    	taul (ms)
    	ghd (mho/cm2)
}

INITIAL {
	rate(v)
	l=linf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ghd = gpeak*l	
	i = ghd*(v-Vrev)			

}


FUNCTION alpt(v(mV)) (1) {			
  alpt = exp(0.0378*zetat*(v-vhtau)) 
}

FUNCTION bett(v(mV)) (1) {			
  bett = exp(0.0378*zetat*gmt*(v-vhtau)) 
}

DERIVATIVE states {     
        rate(v)
        l' =  (linf - l)/taul
}

PROCEDURE rate(v (mV)) { 
        LOCAL a,b,qt
        
        qt=q10^((celsius-temp)/10(degC))	
        a = alpt(v)			
        b = bett(v)
        linf = 1/(1 + exp(-(v-vhakt)/k))	
        taul = b/(qtl*qt*a0t*(1+a))		
}














