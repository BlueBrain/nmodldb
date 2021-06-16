TITLE passive extrasynaptic (i) and synaptic (is) membrane channels 

COMMENT
Both "i" and "is" currents are switched on.
Replicas of NEURON pas.mod membrane mechanism 
for simulation of tonic activation of voltage-independent synaptic 
conductance as excitatory input to passive dendrites.
Reference: Korogod SM and Kulagina IB (1998) Biol Cybern 79:231-240
ENDCOMMENT


UNITS {
   (mV) = (millivolt)
   (mA) = (milliamp)
}

INDEPENDENT { v FROM -100 TO 50 WITH 50   (mV) }

NEURON {
   SUFFIX PasD
   NONSPECIFIC_CURRENT i, is
   RANGE g, erev, gs, es
}

PARAMETER {
   g = 0.0000677254 (mho/cm2)
   erev = -65       (mV)
   gs = 0.0000677254 (mho/cm2)
   es = 0.          (mV)
}

ASSIGNED { i (mA/cm2) is (mA/cm2) }

BREAKPOINT {
   i = g*(v - erev)
   is = gs*(v - es)
}

