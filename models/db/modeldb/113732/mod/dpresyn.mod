: Modified from presyn.mod for dummy lgn cell

COMMENT
calls newspike() if spike occurs
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS DPRESYN
  RANGE in, times, x
}

ASSIGNED {
  times[100]
}

PARAMETER {
  in = 0
  x = 0
}

INCLUDE "presyn.inc"

INITIAL {
  in = 0
  x = 0
}

BREAKPOINT {
  SOLVE pp
}

PROCEDURE pp() {
VERBATIM
  x = fabs(times[(int)in]-t);
ENDVERBATIM
  if (x < 0.05) {
      newspike()
      in = in + 1
    }
}
