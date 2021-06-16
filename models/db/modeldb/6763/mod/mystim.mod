COMMENT
Since this is an electrode current, positive values of i depolarize the cell
and in the presence of the extracellular mechanism there will be a change
in vext since i is not a transmembrane current but a current injected
directly to the inside of the cell.
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS VClamp1
	RANGE vset, gain, i
	ELECTRODE_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	vset = -65  (mV)
	gain = 1 (nA/mV)
	dt  (ms)
        v (mV)
}
ASSIGNED { i (nA) }

BREAKPOINT {
		i = gain*(vset-v)
}
