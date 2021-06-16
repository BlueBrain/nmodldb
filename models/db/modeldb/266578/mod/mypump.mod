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
	SUFFIX mypump
:	USEION Na WRITE iNa VALENCE 1
:        USEION na READ nai
        USEION na WRITE ina
	USEION k  WRITE ik
	RANGE  inapump,ipumpmax,n,km,kk,k,ina,ik,decline,initialdensity,lex,red,lux,green
 
}


PARAMETER {
        dt (ms)
        nai   (mM)
        ipumpmax  = 0.01   (mA/cm2)
        km = 10.0        (mM)
        n=1.5
        kk =  10.0        (mM)
        k = 1.5

        nainit = 4  (mM)
        celsius = 35  (degC)
        T = 1

decline = 0
initialdensity = 0.00208768267   (mA/cm2)
        lex = 100000 (ms)       :   if 100 then we get fast declining pump. if like 100,000. then we dont really get decline        
       lux = 100000 (ms)   
red = 3000 (ms)
  green = 8500 (ms)
:         corrD = 4.79

}

ASSIGNED { 
           ina		(mA/cm2)
           iNa		(mA/cm2)
           ik		(mA/cm2)
        inapump (mA/cm2)
:          ipumpmax (mA/cm2)
        inapumping (mA/cm2)
         xm
         t_wait (ms)
         t_date (ms)
        ipumplux     (mA/cm2)
}


INITIAL{
       ipumpmax = initialdensity
t_wait = 0
}


BREAKPOINT {

if (decline == 1) {

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

ina = 3.0*ipumpmax
ik = -2.0*ipumpmax

}





