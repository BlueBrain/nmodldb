COMMENT

k_slow.mod

voltage gated potassium channel, Hodgkin-Huxley style kinetics.  

Kinetics were fit to data from recordings of nucleated patches derived 
from pyramidal neurons. Data recordings and fits from Alon Korngreen 

I left the voltage shift option untouched. It originally belonged to the 
template of this program. The vshift is set equal to zero. 

Author of the template (na.mod): Zach Mainen, Salk Institute, 1994, 
zach@salk.edu

Author: Sebastian Bergling, MPImF Cell Physiology, 1998,
bergling@mpimf-heidelberg.mpg.de

Author: Alon Korngreen,  MPImF Cell Physiology, 1998,
alon@mpimf-heidelberg.mpg.de

last updated 20/1/1999 by AK

On 5/4/99 AK made some changes to the slow inactivation.  
The changes are not official.  Please use kslow_orig for the official 
version

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX kslow
	USEION k READ ek WRITE ik
	RANGE  a, b, b1,gkslow, gbar,c1,c2
	RANGE  ainf, taua, binf, taub,taub1
	GLOBAL a0, a1, a2, a3, a4, a5, a6
	GLOBAL b0, b11, b2, b3, b4, b5
	GLOBAL bb0,bb1
	GLOBAL v05a, za, v05b, zb
	GLOBAL q10, temp, tadj, vmin, vmax, vshift
}

PARAMETER {
	gbar = 0   	(pS/um2)	: 
	vshift = 0	(mV)		: voltage shift (affects all)
								
	v05a = -15.8	(mV)		: v 1/2 for act (a) 
	za   =  15.5	(mV)		: act slope		
	v05b = -58	(mV)		: v 1/2 for inact (b) 
	zb   = -11  (mV)		: inact slope
		
	a0   =  0.0052  (1/ms*1/mV)		: parameters for alpha and beta for activation
	a1   = 11.1 	(mV)		:      see below
	a2   = 13.1	(mV)		:      see below 
	a3   = 0.01938    (1/ms)		:      see below 
	a4   = -1.27	(mV)		:	see below
	a5   = 71    (mV)
	a6   = -0.0053 (1/ms)	
	
	b0   = 360	(ms)		: fast inact tau (taub) (ms) 
	b11   = 1010	(ms)		:      see below
	b2   = -75	(mV)		:      see below
	b3   = 48	(mV)		:      see below
	b4   = 23.7     (ms/mV)
	b5   = -54      (mV)

	bb0 = 3330	(ms)		: Slow inactivation tau (taub1)
	bb1 = -4	(ms/mV)

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
	gkslow		(pS/um2)
	ek		(mV)
	ainf 		
	binf
	taua (ms)	
	taub (ms)
	taub1 (ms)	
	tadj
}
 

STATE { a b b1 c1 c2 }

INITIAL { 
	trates(v+vshift)
	a = ainf
	b = binf 
	b1= binf
}

BREAKPOINT {
        SOLVE states
        gkslow = tadj*gbar*a*a*(c1*b+c2*b1)
	ik = (1e-4) * gkslow * (v - ek)
} 

LOCAL aexp, bexp,b1exp, z 

PROCEDURE states() {   		:	Computes state variables a and b
        trates(v+vshift) 	:       at the current v and dt.
        a = a + aexp*(ainf-a)
        b = b + bexp*(binf-b)
	b1 = b1 +b1exp*(binf-b1)
        VERBATIM
        return 0;
        ENDVERBATIM
}

PROCEDURE trates(v) { 
                      
    LOCAL tinc

    TABLE ainf, binf, aexp, bexp,b1exp, taua, taub,taub1,c1,c2
	DEPEND dt, celsius, temp, q10, a0, a1, a2, a3, a4, a5, a6, b0, b11, b2, b3, b4,b5, bb0, bb1, v05a, za, v05b, zb
	FROM vmin TO vmax WITH 199

	rates(v): not consistently executed from here if usetable == 1

        tadj = q10^((celsius - temp)/10)
        tinc = -dt * tadj

        aexp = 1 - exp(tinc/taua)
        bexp = 1 - exp(tinc/taub)
	b1exp = 1 - exp(tinc/taub1)
}



PROCEDURE rates(vm) {  

	LOCAL alpha, beta

	alpha=a0*(vm-a1)/(1-exp(-(vm-a1)/a2))
	beta=a3*exp(-(vm-a4)/a5)+a6

	taua=1/(alpha+beta)
	ainf = alpha/(alpha+beta)
	
	taub = b0 + (b11+b4*(vm-b5))*exp(-(vm-b2)*(vm-b2)/(b3*b3))
    	taub1=bb0+bb1*vm
	binf = 1/(1+exp(-(vm-v05b)/zb))
	c1=1-heav(vm+10)*0.55
	c2=heav(vm+10)*0.55
}



FUNCTION heav(y) {
	IF (y<=0) {     
	z=0 
		}
		ELSE {
		z=1
		}
		heav=z 
}
