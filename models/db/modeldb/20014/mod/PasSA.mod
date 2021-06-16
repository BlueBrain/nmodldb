TITLE passive membrane channel

COMMENT
A replica of NEURON pas.mod membrane mechanism 
for passive soma and axon with the parameters as in the models 
described in: 
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
   SUFFIX PasSA
   NONSPECIFIC_CURRENT i
   RANGE g, erev
}

PARAMETER {
   g = .000677254 (mho/cm2)
   erev = -65     (mV)
}

ASSIGNED { i (mA/cm2) }

BREAKPOINT {
   i = g*(v - erev)
}
