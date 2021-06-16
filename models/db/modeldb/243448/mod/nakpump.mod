TITLE Na-K ATPase Pump for bladder small DRG neuron soma model

:Adapted from Tigerholm et al., 2014 and Hamada et al., 2003

: For details refer: 
: A biophysically detailed computational model of bladder small DRG neuron soma 
: Darshan Mandge and Rohit Manchanda, PLOS Computational Biology (2018)

UNITS {
	(mA) = (milliamp)
	(pA) = (picoamp)
	(mV) = (millivolt)
	(S)  = (siemens)
	(molar)  =  (1/liter)
	(mM) =  (millimolar)
	(uF) =   (micro-farad)
	(pF) = (pico-farad)
}

NEURON {
        SUFFIX nakpump
        USEION k READ ko WRITE ik
		USEION na READ nai WRITE ina
        RANGE gbar, ik, ina, ipump, capm
        THREADSAFE
}

PARAMETER {
        gbar = 0.001 (1)	
		capm  (uF/cm2)
		imaxh = 1.62 (pA/pF)
		kh = 6.7 (mM)
		kl = 67.6(mM)
		imaxl = 0.99 (pA/pF)
}

ASSIGNED {
        v	(mV)
        ko  (mM)
        ik  (mA/cm2)
		nai (mM)
		ina (mA/cm2)
		ipump (mA/cm2)
		celsius (degC)
}

BREAKPOINT {        
		ipump = gbar*(imaxh*capm*(1e-3)/(1+(kh/nai)^3) + imaxl*capm*(1e-3)/(1+(kl/nai)^3)) : Hamada et al., 2003

		ina = 3*ipump
		ik = -2*ipump
}
