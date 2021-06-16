NEURON {
  SUFFIX mar
  USEION ca READ cai
  RANGE val, cval
}

ASSIGNED {
  v (millivolt)
  val (millivolt)
	cai (millimolar)
	cval (millimolar)
}

INITIAL {
  val = v
  cval = cai
	
}

BREAKPOINT {
  if (v>val) {
    val = v
  }
  if( cai > cval ) {
	cval = cai
  }
}
