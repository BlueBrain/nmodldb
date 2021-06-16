TITLE PUMP
: FORREST MD (2014) Two Compartment Model of the Cerebellar Purkinje Neuron
 

UNITS {
       (molar) = (1/liter)
        (pA) = (picoamp)
	(mV) =	(millivolt)
        (uS) = (micromho)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}


INDEPENDENT {v FROM -100 TO 50 WITH 50 (mV)}

NEURON {
	SUFFIX newnew
:	USEION Na WRITE iNa VALENCE 1
:        USEION na READ nai
        USEION na READ nai  WRITE ina
	USEION k READ ko  WRITE ik
	RANGE  inapump,ipumpmax,n,km,kk,k,decline,initialdensity,lex,red
 
}


PARAMETER {
        dt (ms)
        nai   (mM)
:        ipumpmax  = 0.04   (mA/cm2)
        km = 0.95        (mM)
        n= 1
        kk =  10.0        (mM)
        k = 1.5

        nainit = 4  (mM)
        celsius = 35  (degC)
        T = 1
decline = 0
initialdensity = 0.0010438413 (mA/cm2)
 lex = 100000 (ms)       :   if 100 then we get fast declining pump. if like 100,000. then we dont really get decline        
   red = 3000 (ms)     
:   corrD = 6.16 
}

ASSIGNED { 
           ina		(mA/cm2)
           iNa		(mA/cm2)
           ik		(mA/cm2)
        inapump (mA/cm2)
          ipumpmax (mA/cm2)
        inapumping (mA/cm2)
         xm
          t_wait (ms)
          ko (mM)
:   corrD
}


INITIAL{
       ipumpmax = initialdensity
t_wait = 0
}


BREAKPOINT {LOCAL fnk

if (decline == 1)
{

if ( t > red) {   
if (t_wait < lex)
{t_wait = t_wait + dt }
else {
ipumpmax = ipumpmax - 0.001
t_wait = 0                    : this resets the counter.
}
}
if (ipumpmax < 0) {ipumpmax = 0}
}


 inapump = ipumpmax*(1/(1 + pow(km/ko,n)))


	ina = 3.0*inapump
	ik = -2.0*inapump
}






