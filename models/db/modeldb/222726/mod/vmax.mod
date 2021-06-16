NEURON {
    SUFFIX vmax
    RANGE vm
}

ASSIGNED {
       v (millivolt)
       vm (millivolt)
}

INITIAL {
    vm = v+70
}

BREAKPOINT { 
   if (v+70>vm) { vm=v+70 }
}
