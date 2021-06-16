: $Id: seed.mod,v 1.1 2005/04/12 05:39:37 billl Exp $
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX nothing
}

: unable to get the drand here to recognize the same fseed used in rand
PROCEDURE vseed (seed) {
  VERBATIM
  srand48((unsigned)_lseed);
  set_seed(_lseed);
  ENDVERBATIM
}

