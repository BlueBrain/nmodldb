TITLE passive membrane channel
                                                                                
UNITS {
        (mV) = (millivolt)
        (mA) = (milliamp)
}
                                                                                
NEURON {
        SUFFIX passsd
        NONSPECIFIC_CURRENT i
        RANGE g, erev, i
}
                                                                                
PARAMETER {
	v		(mV)
        g = .001        (mho/cm2)
        erev = -70      (mV)
}
                                                                                
ASSIGNED { i    (mA/cm2)}
                                                                                
BREAKPOINT {
        i = g*(v - erev)
}

