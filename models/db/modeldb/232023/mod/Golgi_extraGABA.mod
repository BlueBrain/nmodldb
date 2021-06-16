TITLE Extra GABA conductance

COMMENT
Extra GABA conductance in the cerebellar Golgi cell
Written by Sungho Hong and Claus Lang
Computational Neuroscience Unit, Okinawa Institute of Science and Technology, Japan
Supervision: Erik De Schutter

Correspondence: Sungho Hong (shhong@oist.jp)

September 16, 2017
ENDCOMMENT

NEURON {
	POINT_PROCESS Golgi_extraGABA
	NONSPECIFIC_CURRENT i
	RANGE g, i, e

}

PARAMETER {
	g = 216e-06 (microsiemens)
	e = -80 (millivolt)
}
ASSIGNED {
	i (nanoamp)
	v (millivolt)
}
BREAKPOINT {
	i = g*(v - e)
}
