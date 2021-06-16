TITLE Mouse Ventricular Myocyte Model

COMMENT
Corrective temperature-dependent voltage-independent current to maintain resting potential at all temps based up 70 pF cell capacitance
ENDCOMMENT

NEURON {
    SUFFIX ITEMP
    USEION k WRITE ik
    RANGE  ik
}

UNITS {
    (mA) = (milliamp)
}

ASSIGNED {
    celsius (degC)
    ik  (mA/cm2)
}

BREAKPOINT {
     ik = -0.0000322 * (celsius - 43.3(degC))    : (mA/cm2)
}
  