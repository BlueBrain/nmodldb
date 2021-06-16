NEURON {
  SUFFIX max
  RANGE val
}

ASSIGNED {
  v (millivolt)
  val (millivolt)
}

INITIAL {
  val = v
}

AFTER SOLVE {
 if (v > val) {
   val = v
 }
}