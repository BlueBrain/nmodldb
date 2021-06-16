: synaptic current of climbing fiber

NEURON {
    SUFFIX climbing
    RANGE del,tauO,tauC
    NONSPECIFIC_CURRENT  i
    RANGE i, e, g
}

PARAMETER {
    g = 15e-3 (siemens/cm2)  < 0, 1e9 >
    e = 0    (millivolts)
    del = 100 (ms)
    tauO = 0.7 (ms)
    tauC = 6.4 (ms)  

}

ASSIGNED {
    i   (milliamp/cm2)
    v   (millivolt) 
}

INITIAL  { i = 0   }

BREAKPOINT {
    at_time(del)
    if (t < del) {
          i = 0
    } else {
          i = g*(1-exp(-(t-del)/tauO))*exp(-(t-del)/tauC)*(v - e)
    }
}