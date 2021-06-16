COMMENT 

Model for an SODIUM cuRrent recorded in DA neurons.Nav12 style

Activation and inactivation parameters DERIVED FROM MIGLIORE 


ENDCOMMENT


NEURON
 {
	 THREADSAFE
  SUFFIX Na12  USEION na READ ena WRITE ina 
 
RANGE gbar, g, i  

RANGE atau, btau   

RANGE ainf, binf 



 }

UNITS {

(pS) =(picosiemens)
(mV) = (millivolt)
(mA) = (milliamp)

}

PARAMETER { 
  gbar = 50 (pS/microm2)
  Vmid_ac = -28(mV)    : actually corresponds to -20 mV because the m^3 model is rightward shifted
  k_ac = 7.7(mV) 	: corresponds to 7
  k_ina = -10 (mV)
Vmid_ina = -50 (mV)
  
  m=3         
  h=1        : gate parameters according to the HH formalism (m*m*m*h)
celsius = 32 (degC)
q10=1.5

}

 ASSIGNED {
  v	(mV)
  ina 	(mA/cm2)
  i 	(mA/cm2)
  g	(pS/microm2)
  atau (ms)
  btau (ms)
    ainf (1)
  binf (1)
 ena (mV)
 }


STATE {a b}

BREAKPOINT {
  SOLVE states METHOD cnexp

 

  g = (gbar*(a^m)*(b^h)) 
  i = (0.0001)*g*(v-ena)
  ina = i
}

INITIAL {

rates(v)
a= ainf
b=binf


}

DERIVATIVE states {
 rates(v)
  a' = (ainf-a)/atau
  b' = (binf-b)/btau
 
}


FUNCTION a_inf (V (mV)) () {

  a_inf = 1/(1+exp(-(V-Vmid_ac)/k_ac))  : activation system (a*a*a)
}

FUNCTION b_inf (V (mV)) () {
  b_inf = 1/(1+exp(-(V-Vmid_ina)/k_ina)) : inactivation system (b)
}

FUNCTION a_tau (V (mV)) (ms) {
UNITSOFF

a_tau= 0.01+(0.33/(1+((V+20)/30)^2))
: time constant of activation depends on V 

UNITSON
}

FUNCTION b_tau (V (mV)) (ms) {

UNITSOFF

  b_tau = 0.7+(16/(1+((V+50)/8)^2))
 : fast inactivation  


UNITSON 


}



PROCEDURE rates(V (mV)) {
LOCAL qt
UNITSOFF
qt=q10^((celsius-24)/10)
UNITSON 
atau=a_tau(V)/qt

ainf=a_inf(V)

btau=b_tau(V)/qt

binf=b_inf(V)



}
