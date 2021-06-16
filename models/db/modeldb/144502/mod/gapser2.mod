NEURON {
	POINT_PROCESS gap2
	NONSPECIFIC_CURRENT i
	RANGE i, g,t1,t2,g1,g2
	POINTER vgap	
}
PARAMETER {
	v (millivolt)
	g = 0 (nanosiemens) :1 (nanosiemens) : 1nS = 1000 pS
	g1=0 (nanosiemens)
	g2=0 (nanosiemens)
	t1=0 (milliseconds)
	t2=0 (milliseconds)
}
ASSIGNED {
	i (nanoamp)
	vgap (millivolt)
}
BREAKPOINT {
	if(t < t1) {
	     g = g1
	     :printf("1st IF:: t=%g, g=%g\n",t,g)
	}
	if( (t > t1) && (t < t2)) {
	    g = g2
	    :printf("2nd IF:: t=%g, g=%g\n",t,g)
	}
	if(t > t2) {
	    g = g1
	    :printf("3rd IF:: t=%g, g=%g\n",t,g)
	}
	i = (v - vgap)*g*(0.001)
}