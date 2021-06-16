
TITLE Inward-rectifying channel 
 
: Author: Chitaranjan Mahapatra (chitaranjan@iitb.ac.in)
: Computational Neurophysiology Lab
: Indian Institute of Technology Bombay, India 

: For details refer: 
: Mahapatra C, Brain KL, Manchanda R, A biophysically constrained computational model of the action potential 
: of mouse urinary bladder smooth muscle. PLOS One (2018) 





NEURON {
        SUFFIX IH
		 USEION h READ eh WRITE ih VALENCE 1
        RANGE ghbar, gh, vhfa,slp,tf1,tslp,rt,eh,ih
        GLOBAL einf, etau,eh, q10
}

UNITS {
    (mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	
}

PARAMETER {
        v (mV)
        celsius	 (degC)
	    dt (ms)
        ghbar = 0.0001 (mho/cm2) <0,1e9>
        eh = -29 (mV)
		vhfa = -85 (mv)
		slp = 16
		tslp = 10
		rt = 0.5
		q10=2.3
}

:INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}



STATE {
        e
}

ASSIGNED {
	gh (mho/cm2)
	ih (mA/cm2)
	einf
    etau (ms)
}

INITIAL {
    rates(v)
    e = einf
}



BREAKPOINT {
	SOLVE states METHOD cnexp
    
	gh = ghbar*e
    ih = gh*(v - eh)
    }


DERIVATIVE states {     
        rates(v)
        e' = (einf - e)/etau
       	   

}

PROCEDURE rates(v(mV)) { 
             LOCAL a,qt
            qt=q10^((celsius-37)/10)
	
        einf = 1 / (1+exp((v - vhfa) / slp)) 
           
		etau =  rt * qt * (14.617/(1+exp ((49.042+v)/tslp)))
	  


}
FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	       trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}	