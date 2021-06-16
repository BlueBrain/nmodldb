: $Id: precall.mod,v 1.3 1995/12/16 19:00:01 billl Exp $

COMMENT
presynaptic pointer array
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS PRESYN
  RANGE spk_internal, spk, num
  GLOBAL thresh
}

PARAMETER {
  thresh = -20 			: voltage level nec for release
  num = -1	                : number for reporting spk occurrence
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
      VERBATIM
      { extern double hoc_call_func();
      char func[8] = "precall";
      Symbol* s = hoc_lookup(func);
      if (s && num>=0) {
        hoc_pushx(num);
        hoc_call_func(s, 1);
      }
      }
      ENDVERBATIM
    }
  } else { 
    spk_internal = 0            : drop back down at the end of the spike
    spk = 0
  }
}
