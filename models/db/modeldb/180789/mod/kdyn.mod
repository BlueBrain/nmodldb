: simlates()ple first-order model of potassium dynamics
: FORREST MD (2014) Two Compartment Model of the Cerebellar Purkinje Neuron
: from Durstewitz & Gabriel (2006), Cerebral Cortex


NEURON {
        SUFFIX kdyn
        USEION k READ ko,ik WRITE ko 
        RANGE ko, ra, KAF, dep, peak
}

UNITS {
        (mM) = (milli/liter)
        (mA) = (milliamp)
        F    = (faraday) (coul)
}

PARAMETER {
        tck    = 1000   (ms)           : decay time constant
        koinf = 2 (mM)        :3.82 	(mM)      : equilibrium k+ concentration
	kiinf = 140     (mM)	  :
        dep   = 70e-3 (micron)     : depth of shell for k+ diffusion
	KAF   = 0.143 ()		  : K accumulation factor
peak = 3.03 ()
}

ASSIGNED {
        ik     (mA/cm2)
        ra
}

INITIAL {
	ko=koinf
}

STATE { ko (mM) 
}

BREAKPOINT { 
        SOLVE states METHOD derivimplicit
if (ko > peak) { ko = peak}
  if (ko < 2) { ko = 2}
}

DERIVATIVE states {      
 
        ko'= (1e4*(KAF*ik))/(F*dep)     : + (koinf-ko)/tck    
}


