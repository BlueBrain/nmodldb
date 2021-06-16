COMMENT
Leak channel from Manor (Rinzel, Segev, Yarom) 1997
Channel can cause sub-threshold oscillations in interplay with the Ca-current

B. Torben-Nielsen @ HUJI, 7-10-2010
ENDCOMMENT

NEURON {
       SUFFIX leak
       NONSPECIFIC_CURRENT i
       RANGE i,el,gbar
}

UNITS {
      (mV) = (millivolt)
      (mA) = (milliampere)
      (muA) = (microampere)
      (S) = (siemens)
      (mS) = (millisiemens)
}

ASSIGNED {
	 v (mV)
	 i (mA/cm2)
}

PARAMETER {
	  el = -63 (mV)
	  gbar = 0.15 (mS/cm2)
}

BREAKPOINT {
	   i = gbar*(v-el)*(0.001)
}
