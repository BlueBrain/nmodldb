NEURON {
    SUFFIX vmax
    RANGE vm, tm, tstart
}

ASSIGNED {
       v (millivolt)
       vm (millivolt)
       tm (ms)
}

INITIAL {
    vm = v
    tm = 0
}         

PARAMETER{
    tstart = 0   (ms)
    }

BREAKPOINT {
    if (t>tstart){ 
   if (v>vm) { 
   :    if (v > vm + .1){tm =t}
    vm = v
    tm = t
    }
}                 
}
