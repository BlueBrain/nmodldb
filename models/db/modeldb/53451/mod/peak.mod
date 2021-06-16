: peak voltage and time of the peak voltage
NEURON {
	SUFFIX peak
	RANGE tm, vm
}

ASSIGNED {
	tm (ms)
	vm (millivolt)
	v (millivolt)
}

INITIAL {
	tm = t
	vm = v
}

AFTER SOLVE {
	if (v > vm) {
		vm = v
		tm = t
	}
}

