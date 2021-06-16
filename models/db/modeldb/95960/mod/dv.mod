COMMENT
  deltav, deltavmax and vrest
ENDCOMMENT

NEURON {
    THREADSAFE
    SUFFIX dv
    RANGE vrest, deltav, negdeltav, vmax, vmaxt, vmin, vmint
    RANGE dvmax, dvmaxt, dvmin, dvmint, negdvmax
}

ASSIGNED {
	v (millivolt)
	vrest (millivolt)
	deltav (millivolt)
	negdeltav (millivolt)
	vmax (millivolt)
	vmaxt (ms)
	vmin (millivolt)
	vmint (ms)
	dvmax (millivolt)
	dvmaxt (ms)
	dvmin (millivolt)
	dvmint (ms)
	negdvmax (millivolt)
}

INITIAL {
	vrest = v
	deltav = 0
	negdeltav = 0
	vmax = v
	vmaxt = 0
	vmin = v
	vmint = 0
	dvmax = 0
	dvmaxt = 0
	dvmin = 0
	dvmint = 0
	negdvmax = 0
}
BREAKPOINT {
	if (t<dt) {
		vrest=v
	}
	deltav=v-vrest
	negdeltav=vrest-v
	if (v>vmax) {
		vmax=v
		vmaxt=t
	}
	if (v<vmin) {
		vmin=v
		vmint=t
	}
	if (deltav>dvmax) {
		dvmax=deltav
		dvmaxt=t
	}
	if (deltav<dvmin) {
		dvmin=deltav
		dvmint=t
	}
	if (negdeltav>negdvmax) {
		negdvmax=negdeltav
	}
}

