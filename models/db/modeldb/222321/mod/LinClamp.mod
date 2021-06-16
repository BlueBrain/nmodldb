NEURON {
        POINT_PROCESS LinClamp
        RANGE del, t1, tf, amp0, ampf
        ELECTRODE_CURRENT i
}

UNITS {
        (nA) = (nanoamp)
}

PARAMETER {
  del  = 13000   (ms)
  t1   = 33000   (ms)
	tf   = 60000   (ms)
	amp0 = 0.5     (nA)
	ampf = 0.2     (nA)
}

ASSIGNED {
        i (nA)
}

BREAKPOINT {
       if (t < del) {
         i = 0   
}
        else { 
        if (t > del && t < t1) {
          i = ((ampf-amp0)/(t1-del))*(t-del)+amp0
}
        else { 
        if (t > t1 && t < tf) {
          i = ampf
}
}
}
}