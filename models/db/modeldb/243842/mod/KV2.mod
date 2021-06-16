
TITLE Voltage dependent Potassium channel (Kv) 
 
: Author: Chitaranjan Mahapatra (chitaranjan@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: Mahapatra C, Brain KL, Manchanda R, A biophysically constrained computational model of the action potential 
: of mouse urinary bladder smooth muscle. PLOS One (2018) 



NEURON {
	SUFFIX KV2
	USEION k READ ek WRITE ik
	RANGE gkbar,ek, pt,qt,pw,vhi,sli ,vha,sla,ik
    GLOBAL finf,tauf ,ginf,taug, q10
	    
		}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	
    pt = 1
	qt = 1
	pw = 2
	vha = 1.1 (mV)
	sla= 11
	vhi = -58 (mV)
	sli= 15
    gkbar	= 0.006 (mho/cm2)
	ek	= -75	(mV)
	dt          (ms)
	v           (mV)
	celsius = 37(degC)
	q10 = 2.3
	
}

STATE {
	 f g   
}


ASSIGNED {
	
	     ik	   (mA/cm2)
         finf 
		 tauf
		  
		 ginf 
		 taug
		
        
         
}

INITIAL {
      
	  rates(v)
	  f=finf
	  g=ginf
	  	  
       }

BREAKPOINT {
	SOLVE states METHOD cnexp
	        ik = gkbar * f^pw*g* (v - ek)
			
			
  	        }
		
DERIVATIVE states {     
        rates(v)
        g' = (ginf - g)/taug
       f' =  (finf - f)/tauf		   

}

PROCEDURE rates(v(mV)) { LOCAL a,b,qtt
            qtt=q10^((celsius-37)/10)
	
        finf = 1/(1 + exp((vha-v)/sla)) 
        ginf = 1/(1 + exp((v-vhi)/sli))  
	    
          tauf = pt * qtt * (1/(1+(((v+15)/20)^2)))
	       taug = qt * (200+ (10/(1+(((v+54.18)/120)^5))))
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	       trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}	
