TITLE Ornstein-Uhlenbeck process 
: Noisy current with the possibility of adding a bias
: Variance is corrected by dt thus maintaining the second moment of the resulting process
: Written by Patricio Orio

NEURON {
	SUFFIX OU
	RANGE i, D, tau, bias
	NONSPECIFIC_CURRENT i
}

UNITS { (mA) = (milliamp) }

PARAMETER {
	bias = 0 (mA/cm2)
	D = 0.2      (/ms)
	tau = 1 	(ms)
}

ASSIGNED { 
	i (mA/cm2)
	noise (mA/cm2)
	dt (ms)
}

STATE { n (mA/cm2) }

BREAKPOINT {
	SOLVE kin METHOD cnexp
	i = bias + n
}

DERIVATIVE kin {
	noise = 1(mA/cm2) * normrand(0,D*dt)
	n' = (-n + noise)/tau
}
