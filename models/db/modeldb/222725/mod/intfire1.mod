NEURON {
	ARTIFICIAL_CELL IntFire1
	RANGE tau, refrac, m
}

PARAMETER {
	tau = 10 (ms)
	refrac = 5 (ms)
}

ASSIGNED {
	m
	t0(ms)
	refractory
}

INITIAL {
	m = 0
	t0 = t
	refractory = 0 
}

FUNCTION M() {
	if (refractory == 0) {
		M = m*exp(-(t - t0)/tau)
	}else if (refractory == 1) {
		if (t - t0 < .5) {
			M = 2
		}else{
			M = -1
		}
	}
}

NET_RECEIVE (w) {
	if (refractory == 0) {
		m = m*exp(-(t - t0)/tau)
		t0 = t
		m = m + w
		if (m > 1) {
			refractory = 1
			m = 2
			net_send(refrac, refractory)
			net_event(t)
		}
	}else if (flag == 1) {
		t0 = t
		refractory = 0
		m = 0
	}
}
