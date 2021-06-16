TITLE Leakm current
:NOTE 1S=1mho Neuron wants the units in mhos not millisiemens, please not the conversion! 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX Leakm
        NONSPECIFIC_CURRENT il
        RANGE  gl, el
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
       
        gl = .00001 (mho/cm2) : changed from 0.075 mS by KMB
        el = -75 (mV)
}
  
ASSIGNED {
	 v (mV)
        il (mA/cm2)
}
 
BREAKPOINT {
:        SOLVE states
        il = gl*(v - el)
}
 
