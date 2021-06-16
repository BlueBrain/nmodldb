NEURON {

    POINT_PROCESS Gap
    RANGE gmax, i 
    POINTER vgap
    NONSPECIFIC_CURRENT i
}

UNITS {

  (nA) = (nanoamp)
  (mV) = (millivolt)
  (nS) = (nanosiemens)
}

PARAMETER { gmax = 0 (nanosiemens) }
    
ASSIGNED {

    v    (mV)
    vgap (mV)
    i    (pA)
}
 
BREAKPOINT { 

  if (gmax>0) {i = gmax * (v-vgap) }

}
