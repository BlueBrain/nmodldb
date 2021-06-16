COMMENT 

Model for an kV4 cuRrent recorded in DA neurons.
This current has two inactivation rates and a rapid rate for recovery from inactivation
Activation and inactivation parameters are from amendola

ENDCOMMENT


NEURON
 {
	 THREADSAFE
  SUFFIX kaDasoma  USEION k READ ek WRITE ik 
 
RANGE gbar, g, i  

RANGE atau, btau

RANGE ainf, binf

RANGE taurecov


 }

UNITS {

(pS) =(picosiemens)
(mV) = (millivolt)
(mA) = (milliamp)

}

PARAMETER { 
  gbar = 50 (pS/microm2)
  Vmid_ac =-30 (mV)
  k_ac = 7 (mV) 
  Vmid_ina = -75 (mV) 
  k_ina = -7 (mV)
   taurecov=25 (ms)
  Vshift=-90 (mV) : potential at which inactivation rates are replaced by taurecov
   m=1         
  h=1        : gate parameters according to the HH formalism (m*m*m*h)
mean_inac = 32 (ms) :  inactivation duration 
}

 ASSIGNED {
  v	(mV)
  ik 	(mA/cm2)
  i 	(mA/cm2)
  g	(pS/microm2)
  atau (ms)
  btau (ms)
   
  ainf (1)
  binf (1)
 
ek (mV)
 }


STATE {a b b2}

BREAKPOINT {
  SOLVE states METHOD cnexp

 
  g = gbar*(a^m)*(b^h)
  i = (0.0001)*g*(v-ek)
  ik = i
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

a_tau= 1.029 + (4.83/(1+exp((V+57)/6.22)))
: time constant of activation depends on V 

UNITSON
}

FUNCTION b_tau (V (mV)) (ms) {

UNITSOFF

  b_tau = taurecov + (mean_inac+(78.4/(1+exp(V+68.5)/5.95))-taurecov)/(1+exp((-V+Vshift)*5))
 : fast inactivation  


UNITSON 

}



PROCEDURE rates(V (mV)) {





atau=a_tau(V)

ainf=a_inf(V)

btau=b_tau(V)

binf=b_inf(V)



}
