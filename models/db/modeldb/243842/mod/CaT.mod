
TITLE T-type calcium  channel  

: Author: Chitaranjan Mahapatra (chitaranjan@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: Mahapatra C, Brain KL, Manchanda R, A biophysically constrained computational model of the action potential 
: of mouse urinary bladder smooth muscle. PLOS One (2018) 





NEURON {
	SUFFIX CaT
	
	USEION ca  READ eca WRITE ica VALENCE 2
	
         
	    RANGE gcat, eca, mt,ht,x,y,vhfa,vhfi,sla,sli,icat
        GLOBAL cinf,binf,tauc,taub, q10
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}


PARAMETER {
	  v (mV)
	  celsius	(degC)
	  dt (ms)
	  gcat= 0.0002(mho/cm2)
       mt = 0.1 
	   ht = 0.5
	   x = 2
	   y = 1
	   	   
	   vhfa= -36.9 (mV) 
	   vhfi = -74.8 (mV) 
	   sla = 6.6 (mV) 
	   sli = 50 (mV) 
	     
		 q10=2
		
		
	eca = 51 (mV)
}

ASSIGNED {
	   
		ica 		(mA/cm2)
		icat 		(mA/cm2)
        cinf
        binf      
        tauc
        taub
		
}


STATE {
	    b
        c
}


INITIAL {
	rates(v)
	b=binf
	c=cinf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	    
	ica = gcat*b^x*c^y*(v-eca)
	icat = ica
	
	}


DERIVATIVE states { 
        rates(v)
        b' = (binf - b)/taub
        c' =  (cinf - c)/tauc
}

PROCEDURE rates(v (mV)) { 

        LOCAL a,qt
		UNITSOFF
        qt=q10^((celsius-37)/10)
			
						 
	 binf = 1/(1 + exp((vhfa-v)/sla)) 
	 cinf = 1/(1 + exp((v-vhfi)/sli)) 
		
	 taub = mt* qt *(0.45 + (3.9/(1+((v +66)/26)^2))) 
	 tauc = ht* qt * (150 - (150/((1+ exp((v-417.43)/203.18))*(1+exp(-(v+61.11)/8.07))))) 
	
	UNITSON
	          
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}	