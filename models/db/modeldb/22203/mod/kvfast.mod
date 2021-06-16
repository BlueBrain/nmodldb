COMMENT

modified from 
k_fast.mod  with the purpose to "symbolify"...   A.S.Nov 99

voltage gated potassium channel, Hodgkin-Huxley style kinetics.  

Kinetics were fit to data from recordings of nucleated patches derived 
from pyramidal neurons. Data recordings and fits from Alon Korngreen 

Inactivation time constant is 4 times faster then data.

I left the voltage shift option untouched. It originally belonged to the 
template of this program. The vshift is set equal to zero. 

Author of the template (na.mod): Zach Mainen, Salk Institute, 1994, 
zach@salk.edu

Author: Alon Korngreen, MPImF Cell Physiology, 1998,
bergling@mpimf-heidelberg.mpg.de

Last modified by AK 17.1.1999

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kvfast
	USEION k READ ek WRITE ik
	RANGE  a, b, gkvfast, gbar
	RANGE  ainf, taua, binf, taub
:	GLOBAL ainf, taua, binf, taub
	GLOBAL a0, a1, a2, a3, a4, a5, a6
	GLOBAL b0, b1, b3, b4, b5, b6
	GLOBAL q10, temp, tadj, vmin, vmax, vshift
}

PARAMETER {
	gbar = 0   	(pS/um2)	: 
	vshift = 0	(mV)		: voltage shift (affects all)
									
	a0   =  3.4  	(1/ms)		: activation alpha   
	a1   =  10 	(mV)		:      see below
	a2   =  30 	(mV)		:      see below 

	a3   = 0.7 	(1/ms)		: activation beta
	a4   = 6    	(mV)		:   
	a5   = 180  	(mV)		:      see below
	a6   = -0.45 	(1/ms)		:      see below 

	b0   = 0.0001          (1/ms)		:  inactivation alpha
	b1   = 0.06	(1/mV)		:      see below
	


	b3   = 0.12    (1/ms)		:   inactivation beta
	b4   = -50      (mV)
	b5   = 8         (mV)
	b6   = 0.002  (1/ms)
	

	temp = 21	(degC)		: original temp 
	q10  = 2.3			: temperature sensitivity

	v 		(mV)
	dt		(ms)
	celsius		(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ik 		(mA/cm2)
	gkvfast		(pS/um2)
	ek		(mV)
	ainf 		
	binf
	taua (ms)	
	taub (ms)	
	tadj
}
 

STATE { a b }

INITIAL { 
	trates(v+vshift)
	a = ainf
	b = binf 
}

BREAKPOINT {
        SOLVE states
        gkvfast = tadj*gbar*a*a*a*a*b
	ik = (1e-4) * gkvfast * (v - ek)
} 

LOCAL aexp, bexp, z 

PROCEDURE states() {   		:	Computes state variables a and b
        trates(v+vshift) 	:       at the current v and dt.
        a = a + aexp*(ainf-a)
        b = b + bexp*(binf-b)
        VERBATIM
        return 0;
        ENDVERBATIM
}

PROCEDURE trates(v) { 
                      
    LOCAL tinc

    TABLE ainf, binf, aexp, bexp, taua, taub
	DEPEND dt, celsius, temp, q10, a0, a1, a2, a3, a4, a5,a6, b0, b1, b3, b4, b5, b6
	FROM vmin TO vmax WITH 199

	rates(v): not consistently executed from here if usetable == 1

        tadj = 1 :  q10^((celsius - temp)/10)
        tinc = -dt * tadj

        aexp = 1 - exp(tinc/taua)
        bexp = 1 - exp(tinc/taub)
}



PROCEDURE rates(vm) {  

	LOCAL alphaA, betaA,AlphaI,BetaI
	
	alphaA=a0/(1+exp(-(vm-a1)/a2))
	betaA=a3*exp(-(vm-a4)/a5)+a6
    	taua = 1/(alphaA+betaA)
	ainf = alphaA/(alphaA+betaA)


	AlphaI=b0*exp(-b1*vm)
	BetaI=b3/(1+exp(-(vm-b4)/b5))+b6
    	taub = 0.25/(AlphaI+BetaI)
	binf = AlphaI/(AlphaI+BetaI)
}


