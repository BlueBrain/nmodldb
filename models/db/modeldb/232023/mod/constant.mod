TITLE Constant current

COMMENT
Constant current
Written by Sungho Hong and Claus Lang
Computational Neuroscience Unit, Okinawa Institute of Science and Technology, Japan
Supervision: Erik De Schutter

Correspondence: Sungho Hong (shhong@oist.jp)

September 16, 2017
ENDCOMMENT

NEURON {
    SUFFIX constant
    NONSPECIFIC_CURRENT i
    RANGE i,ic
}

UNITS {
    (mA) = (milliamp)
}

PARAMETER {
    ic = 0 (mA/cm2)
}

ASSIGNED {
    i (mA/cm2)
}

INITIAL {
    i=ic
}

BREAKPOINT {
    i=ic
}

