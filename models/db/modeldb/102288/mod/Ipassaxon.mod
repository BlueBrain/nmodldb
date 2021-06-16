TITLE passive membrane channel
                                                                                
UNITS {
        (mV) = (millivolt)
        (mA) = (milliamp)
}
                                                                                
NEURON {
        SUFFIX passaxon
        NONSPECIFIC_CURRENT i
        RANGE g, erev
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

