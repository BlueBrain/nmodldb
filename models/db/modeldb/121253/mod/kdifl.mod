COMMENT
Longitudinal diffusion of potassium
(equivalent modified euler with standard method and
equivalent to diagonalized linear solver with CVODE )
ENDCOMMENT

NEURON {
	SUFFIX kdifl
	USEION k READ ik WRITE ko
	RANGE ko, Dk, sp, rseg
}

UNITS {
  	(molar) = (1/liter)
  	(mM) = (millimolar)
	(um) = (micron)
	(mA) = (milliamp)
	FARADAY = (faraday) (coulomb)
	PI = (pi) (1)
}

PARAMETER {
	Dk = 1.85 (um2/ms)
	sp = 0.01 (um)
	rseg = 1.7 (um)
	extracellularVolumePerLength (um2)
	crossSectionalArea (um2)

}

ASSIGNED { ik (mA/cm2) }

STATE { ko (mM) }

INITIAL {
	extracellularVolumePerLength = PI * ( (rseg+sp)^2 - rseg^2 )
	crossSectionalArea = extracellularVolumePerLength
	
}

BREAKPOINT { SOLVE conc METHOD sparse }

KINETIC conc {
	COMPARTMENT extracellularVolumePerLength {ko}
	LONGITUDINAL_DIFFUSION Dk*crossSectionalArea {ko}
	~ ko << (ik/(FARADAY)*2*PI*rseg*(1e4))
}
