TITLE lk.mod   leak channels
 
COMMENT
This is the original Hodgkin-Huxley treatment for the set of leak channels found
in the squid giant axon membrane.
Some parameters have been changed to correspond to McIntyre and Grill (2002) "Extracellular
stimulation of central neurons"
Author: Balbi
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
		(S) = (siemens)
}
 
NEURON {
        SUFFIX lk
        NONSPECIFIC_CURRENT il
        RANGE gl, el
}
 
PARAMETER { : values for IS
        gl = 0.007 (S/cm2)	  <0,1e9>
        el = -70 (mV)
}
 
ASSIGNED {
        v (mV)

        il (mA/cm2)
}
 
BREAKPOINT {
        il = gl*(v - el)
}
