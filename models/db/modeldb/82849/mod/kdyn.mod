: simple first-order model of potassium dynamics

NEURON {
        SUFFIX kdyn
        USEION k READ ko,ik WRITE ko 
        RANGE ko, ki
	
}

UNITS {
        (mM) = (milli/liter)
        (mA) = (milliamp)
        F    = (faraday) (coul)
}

PARAMETER {
        tk    = 7   (ms)           : decay time constant
        koinf = 3.82 	(mM)      : equilibrium k+ concentration
	kiinf = 140     (mM)	  :
        dep   = 70e-3 (micron)     : depth of shell for k+ diffusion
	KAF   = 2 ()		  : K accumulation factor

	 
}

ASSIGNED {
        ik     (mA/cm2)
}


STATE { ko (mM)  }

BREAKPOINT { 
        SOLVE states METHOD cnexp
}

DERIVATIVE states {
    
        ko'= 1e4*KAF*ik/(F*dep) + (koinf-ko)/tk
 
}





