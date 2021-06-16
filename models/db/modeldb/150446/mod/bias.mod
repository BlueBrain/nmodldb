TITLE Current bias used as alternative to IClamp to bias voltage of many neurons
    NEURON {
      SUFFIX bias
      NONSPECIFIC_CURRENT i
      RANGE i, amp
    }
    UNITS {
      (mA) = (milliamp)
    }
    PARAMETER {
      amp = 0 (mA/cm2)
    }
    ASSIGNED {
      i (mA/cm2)
    }
    BREAKPOINT {
      i = amp
    }