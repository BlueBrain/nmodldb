NEURON {
	SUFFIX nothing
}

PROCEDURE flushf() {
  VERBATIM
    fflush(NULL);
  ENDVERBATIM
}
