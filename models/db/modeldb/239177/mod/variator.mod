NEURON {
	THREADSAFE
	POINT_PROCESS variator
	RANGE t0, t1, a0, a1, a
	POINTER var
}
PARAMETER {
	t0 ()
	t1 ()
	a0 ()
	a1 ()
}
ASSIGNED {
	a ()
	var ()
}
BREAKPOINT {
	if( t < t0 ){
		a = a0
	}else{
		if( t > t1 ){
			a = a1
		}else{
			a = a0 + (a1-a0)*(t-t0)/(t1-t0)
		}
	}
	var = a
}

