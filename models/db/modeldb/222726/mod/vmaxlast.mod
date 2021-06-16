NEURON {
    SUFFIX vmaxlast
    RANGE vm
}

ASSIGNED {
       v (millivolt)
       vm (millivolt)
}

PARAMETER {
	tcheck=20710 (ms)
}

INITIAL {
    vm = v+70
}

BREAKPOINT { 
	if (t>tcheck) {
	   if (v+70>vm) {vm=v+70 }
	}
}
