COMMENT 

Model for an kVdr current recorded in DA neurons.
This current has a form of n4.No inactivation
Activation  parameters are derived from Schild et al  J physiol (1993)  see also paton 

ENDCOMMENT


NEURON
 {
	 THREADSAFE
  SUFFIX kdrDA  USEION k READ ek WRITE ik 
 
RANGE gbar, g, i  

RANGE atau

RANGE ainf 


 }

UNITS {

(pS) =(picosiemens)
(mV) = (millivolt)
(mA) = (milliamp)

}

PARAMETER { 
  gbar = 20 (pS/microm2)
  ek = -90 (mV)
  Vmid_ac = -30(mV)
  k_ac = 9 (mV)     
n=4         
         : gate parameters according to the HH formalism (n*n*n*n)
celsius = 32 (degC)
tau_act=4 (ms)
q10=1.5
}

 ASSIGNED {
  v	(mV)
   ik 	(mA/cm2)
  i 	(mA/cm2)
  g	(pS/microm2)
  atau (ms)
  ainf (1)
  
 }


STATE {a}

BREAKPOINT {
  SOLVE states METHOD cnexp

 
  g = gbar*(a^n)
  i = (0.0001)*g*(v-ek)
  ik = i
}

INITIAL {

rates(v)
a= ainf

}

DERIVATIVE states {
 rates(v)
  a' = (ainf-a)/atau
  
}


FUNCTION a_inf (V (mV)) () {

  a_inf = 1/(1+exp(-(V-Vmid_ac)/k_ac))  

: activation system (a*a)
}


FUNCTION a_tau (V (mV)) (ms) {
UNITSOFF

a_tau= (tau_act * exp(-(0.000729)*((V +32)^2))) + 4 

: time constant of activation depends on V 

UNITSON
}



PROCEDURE rates(V (mV)) {
LOCAL qt

qt=q10^((celsius-24)/10)

atau=a_tau(V)/qt

ainf=a_inf(V)


}
