: $Id: presyn.mod,v 1.4 1995/10/21 23:24:50 billl Exp $

COMMENT
presynaptic pointer array
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS PRESYN
  RANGE spk_internal, spk
  GLOBAL thresh
}

PARAMETER {
  thresh = 0 			: voltage level nec for release
}

ASSIGNED { 
  spk                           : available for user monitoring of spiking
  spk_internal                  : internal use only (if taken externally ...)
  v
}

INCLUDE "presyn.inc"

INITIAL {
  spk = 0
  spk_internal = 0
}

BREAKPOINT {
  SOLVE pp
}

PROCEDURE pp() {
  if (v > thresh) {    
    if (spk_internal == 0) {  
      newspike()                : only allow this to happen once
      spk_internal = 1
      spk = 1
    }
  } else { 
    spk_internal = 0            : drop back down at the end of the spike
    spk = 0
  }
}
