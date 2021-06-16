COMMENT
 A gap junction connecting two compartments.
 1: Create one gapJ object in each of the compartments which 
    are to be connected.
 2: Connect the two with the command
    setpointer gJ1.vgap, c2.comp2.v( 0.5 )
    setpointer gJ2.vgap, c1.comp1.v( 0.5 )
    where gJ1 is the gapJ object in compartment comp1 in cell c1, and
    gJ2 is the gapJ object in compartment comp2 in cell c2
 3: Set gJ1.g = gJ2.g = g0

 Author: Fredrik Edin, 2003
 Address: freedin@nada.kth.se

ENDCOMMENT

NEURON {
	POINT_PROCESS gapJ
	NONSPECIFIC_CURRENT i
	RANGE g, i
	POINTER vgap
}

PARAMETER {
	g	(umho)
}

ASSIGNED {
	i 	(milliamp)
	v 	(millivolt)
	vgap 	(millivolt)
}

BREAKPOINT {
	i = g*(v - vgap)
}

