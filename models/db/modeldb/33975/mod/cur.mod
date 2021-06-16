TITLE steady current, amp is positive inward

UNITS {
	(mA) = (milliamp)
}

NEURON {
	SUFFIX cur
	NONSPECIFIC_CURRENT i
	RANGE amp
}

PARAMETER {
	amp (milliamp/cm2)
}

ASSIGNED { i (mA/cm2)}

BREAKPOINT {
	i = -amp
}
