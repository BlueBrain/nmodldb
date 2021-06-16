NEURON {
    SUFFIX gScale
    RANGE gVal,vX1,vX2,vCP,dX2,dCP
}

ASSIGNED {
	gVal : scaled conductance value
	vX1 : voltage with one synapse
	vX2 : voltage with double synapse
	vCP : voltage with compound synapse
	dX2 : deviation from linearity, double synapse
	dCP : deviation from linearity, compound synapse
}