COMMENT
Longitudinal and radial diffusion of potassium
ENDCOMMENT

NEURON {
	SUFFIX kdifrl
	USEION k READ ik WRITE ko
	RANGE ko, Dk, sp
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
	nodalGap = 1.9 (um)
	rn = 0.95 (um)
	extracellularVolumePerLength (um2)
	crossSectionalArea (um2)
	kbathVolumePerLength = 1e15 (um2)
	radialCrosssectionalAreaPerLength (um)
	radialRateConstant (um2/ms)
	kbath = 3 (mM)
}

ASSIGNED { ik (mA/cm2) }

STATE { ko (mM) }

INITIAL {
	extracellularVolumePerLength =PI * ( (rn+nodalGap)^2 - rn^2 )
	crossSectionalArea = extracellularVolumePerLength

	radialCrosssectionalAreaPerLength = 2 * PI * ( rn + nodalGap )
	radialRateConstant = Dk * radialCrosssectionalAreaPerLength / nodalGap
}

BREAKPOINT { SOLVE conc METHOD sparse }

KINETIC conc {
	COMPARTMENT crossSectionalArea { ko }
	COMPARTMENT kbathVolumePerLength { kbath }
	LONGITUDINAL_DIFFUSION Dk*extracellularVolumePerLength { ko }
	~ ko << (ik/(FARADAY)*2*PI*rn*(1e4))
	~ ko <-> kbath (radialRateConstant,radialRateConstant)
}
