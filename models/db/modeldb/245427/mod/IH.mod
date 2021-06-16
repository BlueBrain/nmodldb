COMMENT 

Model for an Ih amendola
ENDCOMMENT


NEURON
 {
	 THREADSAFE
  SUFFIX Ih
NONSPECIFIC_CURRENT i
 
RANGE gbar, g, i  

RANGE tauh
RANGE ainf


 }

UNITS {

(pS) =(picosiemens)
(mV) = (millivolt)
(mA) = (milliamp)

}

PARAMETER { 
  gbar = 2(pS/microm2)
  eh = -40 (mV)
  Vmid_ac = -92 (mV)
  k_ac = -7.25 (mV) 
  m=1         
  h=0        : gate parameters according to the HH formalism (m*m*m*h)
celsius = 32 (degC)
q10=1.5
}

 ASSIGNED {
  v	(mV)
   i 	(mA/cm2)
  g	(pS/microm2)
 tauh (ms)
   ainf (1)
binf(1)
 
 }


STATE {a b}

BREAKPOINT {
  SOLVE states METHOD cnexp
  g = (gbar*(a^m)*(b^h)) 
  i = (0.0001)*g*(v-eh)
  
}

INITIAL {

rates(v)
a= ainf
b=binf


}

DERIVATIVE states {
 rates(v)
  a' = (ainf-a)/tauh
 
  
}


FUNCTION a_inf (V (mV)) () {

  a_inf = 1/(1+exp(-(V-Vmid_ac)/k_ac))  : activation system 
}



FUNCTION a_tauh (V (mV)) (ms) {
UNITSOFF

a_tauh= 556+ 1100*exp(-0.5*((V)/11.06)^2)
: time constant of activation depends on V 

UNITSON
}




PROCEDURE rates(V (mV)) {

LOCAL qt

qt=q10^((celsius-24)/10)

tauh=a_tauh(V)/qt

ainf=a_inf(V)


}
