: Ksub sub-threshold K current with steep voltage activation

NEURON {
     SUFFIX ksub
     USEION k READ ek WRITE ik
     RANGE gk, ik
}

UNITS {
     (S)  = (siemens)
     (mV) = (millivolt)
     (mA) = (milliamp)
}

PARAMETER { gk = 3e-5 (S/cm2) }

ASSIGNED {
     v       (mV)
     ek      (mV)
     ik      (mA/cm2)
}


BREAKPOINT {
     ik = gk * ninf(v)^3 * (v - ek)
}


FUNCTION ninf (Vm (mV)) {
     UNITSOFF
     ninf = 1 / (1 + exp(-(Vm+44.5)/3)) 
     UNITSON
}



