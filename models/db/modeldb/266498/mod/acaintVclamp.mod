: Intracellular calcium ion accumulation
NEURON {
  SUFFIX acaint
  USEION ca READ ica WRITE cai
  GLOBAL nb
  RANGE a, ku, kr, Bi, diffOc
}
UNITS {
  (mV)    = (millivolt)
  (mA)    = (milliamp)
  FARADAY = (faraday) (coulombs)
  (molar) = (1/liter)
  (mM)    = (millimolar)
}
PARAMETER {
  a = 15e-4 (cm)  : radius of cell
  ku = 100 (/mM/ms) : rate constant for calcium buffer binding
  kr = 0.238 (/ms) : rate constant for calcium buffer release
  nb = 4 : number of binding sites on Calmodulin
  Bi = 0.001 (mM) :Concentration of Calmodulin
  SA = 2.82743E-05 (cm2) :Surface area of cell
  Vol = 1.27e-8 (cm3) :vol of cell
}
ASSIGNED { 

ica  (mA/cm2)
diffOc (/ms)


 }
STATE { cai  (mM) Oc } :Oc is the fraction of Calmodulin binding sites that are occupied
BREAKPOINT { SOLVE state METHOD derivimplicit }
INITIAL {

Oc=.05

}
DERIVATIVE state {
  LOCAL diffOc
  diffOc=ku*cai*(1-Oc)-kr*Oc :Differential equation governing Calmodulin binding
  Oc'=diffOc
  cai' = -ica*(SA)/Vol/(2*FARADAY) + -nb*Bi*diffOc  
}