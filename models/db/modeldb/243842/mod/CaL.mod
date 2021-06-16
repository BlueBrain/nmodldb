
TITLE L-type calcium  channel  

: Author: Chitaranjan Mahapatra (chitaranjan@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: Mahapatra C, Brain KL, Manchanda R, A biophysically constrained computational model of the action potential 
: of mouse urinary bladder smooth muscle. PLOS One (2018) 




NEURON {
	SUFFIX CaL
	USEION ca READ eca WRITE ica VALENCE 2
	
	RANGE  gcal, et,ct,mt,x, y,ki, ical:,cai,cao,ical
	RANGE vhfa,slpa,vhfi,slpi
	GLOBAL cinf,ctau,cvinf,cvtau,dcinf,dctau,q10, eca
	
	
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
         FARADAY = (faraday) (coulomb)
	 R = (k-mole) (joule/degC)
} 

PARAMETER {
	
	 gcal = 0.0004   	(mho/cm2)
	 celsius		(degC)
	 v 		(mV)
	 ki	= .001 	(mM)
     et = 1
	 ct = 1
	 mt = 1
     vhfa= -9   (mV)
	 slpa = 8    (mV)
	 vhfi = 30 (mV)
	 slpi = 13    (mV)
	        
		   
	 q10=2.3
	 cai = 0.00015 (mM)
     
	eca = 51 (mV)
	x = 1
	y = 4
		   		   
}


ASSIGNED {
	ica 		(mA/cm2)
	ical 		(mA/cm2)
	cinf
	ctau (ms)
	cvinf
	cvtau (ms) 
	dctau (ms)
	dcinf
	:eca (mV)
	
}
 

STATE { 
        c 
        cv 
		dc 
	}


INITIAL {
	rates(v)
	c=cinf  
	cv=cvinf
	dc = dcinf
	
}

BREAKPOINT {
        SOLVE states METHOD cnexp
       
	   ica = (gcal*c^x*cv*dc*(v-eca)) 
	                         
		ical = ica			
							 
} 



DERIVATIVE states {   
         rates(v)      
          c' = (cinf-c)/ctau
          cv' = (cvinf-cv)/cvtau
          dc' = (dcinf-dc)/dctau

}

PROCEDURE rates(v (mV)) { :callable from hoc
      LOCAL a,qt
      qt=q10^((celsius-37)/10)         
	         
	  cinf = 1/ (1 + exp ((vhfa-v)/slpa)) 
	  	     
	 ctau = mt* qt *(0.001/(1+exp (-(v+22)/30)))	

          
	  cvinf = 1/(1+exp((vhfi+v)/slpi)) 
		     
	  dcinf = ki / (ki + cai)
	  dctau = ct * 20 
             		 
	  cvtau =  et* (90 *(1.0-(1.0/(1.0+exp((v+14)/45)))*(1.0/(1.0+exp(-(v+9.8)/8.89)))))	 		 
		 			 
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	       trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}	