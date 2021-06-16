NEURON {
	POINT_PROCESS stam
	RANGE d
}

ASSIGNED {
	d
}

FUNCTION draw(lambda) {
	d = exprand(lambda)
}
