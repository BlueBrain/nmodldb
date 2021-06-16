NEURON	{ 
  POINT_PROCESS scopRandom
  RANGE r
}

PARAMETER {
	s
}

ASSIGNED {
	r
}

PROCEDURE seed(x) {
	s = x
	set_seed(s)
}

PROCEDURE next() {
	r = scop_random()
}

INITIAL {
	r = 0
}	

