TITLE passive extrasynaptic (i) and synaptic (is) membrane channels

COMMENT
Current "i" is off, current "is" is switched on.
Replicas of NEURON pas.mod membrane mechanism 
for simulation of tonic activation of voltage-independent synaptic 
conductance as excitatory input to active dendrites.
References: 
1. Korogod SM and Kulagina IB (1998) Biol Cybern 79:231-240
2. Korogod SM, Kulagina IB, and Tyc-Dumont S (1998) Neirofiziologiya/Neurophysiology, 
   Vol.30, Nos.4/5, pp.259-264
   (Kluwer Academic/ Plenum Publishers English version: 
   Neurophysiology 30(4.5):203-207, 1999)
ENDCOMMENT

UNITS {
   (mV) = (millivolt)
   (mA) = (milliamp)
}

INDEPENDENT { v FROM -100 TO 50 WITH 50   (mV) }

NEURON {
   SUFFIX PasS
   NONSPECIFIC_CURRENT i, is
   RANGE g, erev, gs, es
}

PARAMETER {
   g = 0.0           (mho/cm2)
   erev = -65        (mV)
   gs = 0.0000677254 (mho/cm2)
   es = 0.0          (mV)
}

ASSIGNED { i (mA/cm2) is (mA/cm2) }

BREAKPOINT {
   i = g*(v - erev)
   is = gs*(v - es)
}
