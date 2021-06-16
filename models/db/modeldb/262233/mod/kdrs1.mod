COMMENT
Conceptual model:  Delayed rectifier current for 
  a model of a fast-spiking cortical interneuron.

Authors and citation:
  Golomb D, Donner K, Shacham L, Shlosberg D, Amitai Y, Hansel D (2007).
  Mechanisms of Firing Patterns in Fast-Spiking Cortical Interneurons. 
  PLoS Comput Biol 3:e156.

Original implementation and programming language/simulation environment:
  by Golomb et al. for XPP
  Available from http://senselab.med.yale.edu/modeldb/ShowModel.asp?model=97747

This implementation is by N.T. Carnevale and V. Yamini for NEURON.
ENDCOMMENT

NEURON {
  SUFFIX kdrs1
  USEION k READ ek WRITE ik
  RANGE g,gkdr,ikk
}

UNITS {
  (S) = (siemens)
  (mV) = (millivolt)
  (mA) = (milliamp)
}

PARAMETER {
  gkdr = 0.225 (S/cm2)
  theta_hn = -20 (mV):-20
  sigma_n = 10.4 (mV) : 6.8
  q10=3
}

ASSIGNED {
  v (mV)
  ek (mV)
  ik (mA/cm2)
  g (S/cm2)
  qt
  ikk (mA/cm2)
}

STATE {n}

BREAKPOINT {
  SOLVE states METHOD cnexp
  g = gkdr * n*n
  ik = g * (v-ek)
  ikk= ik
}

INITIAL {
  n = ninfi(v)
}

DERIVATIVE states {
  n' = (ninfi(v)-n)/taun(v)
}

FUNCTION ninfi(v (mV)) {
  UNITSOFF
  ninfi=1/(1 + exp(-(v-theta_hn)/sigma_n))
  UNITSON
}

FUNCTION taun(v (mV)) (ms) {
  UNITSOFF
  :qt = q10^((celsius - 24)/10)
  qt = q10^((celsius - 37)/10)

  
 : taun = (((0.087 + 85.4 / (1 + exp ((v+14.6)/8.6))) * (0.087 + 55.4 / (1 + exp (-(v-1.3)/18.7)))))  
 : taun = (((0.087 + 11.4 / (1 + exp ((v+12.6)/18.6))) * (0.087 + 25.4 / (1 + exp (-(v-18.3)/18.7))))) 
  taun = (((0.087 + 17.4 / (1 + exp ((v+35.6)/9.6))) * (0.087 + 25.4 / (1 + exp (-(v-1.3)/18.7))))) : equation from kv2 paper


  
  taun=taun/qt
  UNITSON
}
