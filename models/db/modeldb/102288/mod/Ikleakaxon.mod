TITLE potassium leak channel
                                                                                
UNITS {
        (mV) = (millivolt)
        (mA) = (milliamp)
}
                                                                                
NEURON {
        SUFFIX Kleakaxon
        USEION k READ ek WRITE ik
        RANGE gkl
}
                                                                                
PARAMETER {
	v		(mV)
        gkl = .001        (mho/cm2)
        ek = -70      (mV)
}
                                                                                
ASSIGNED { ik    (mA/cm2)}
                                                                                
BREAKPOINT {
        ik = gkl*(v - ek)
}

