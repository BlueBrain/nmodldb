: 3-D coordinates for each segment 

NEURON {
	   SUFFIX co
	   RANGE x, y, z, pt0, vdist
}

ASSIGNED {
		 x(micron)
		 y(micron)
		 z(micron)
		 pt0 :number of the 3d point just before the segment midpoint
		 vdist(micron) :projection of segment midpoint distance to an origin on a vortex 
}
