NEURON {
	SUFFIX xtra
	RANGE rx, is
	RANGE x, y, z, XX, YY, DX,DY, FX, FY, Exs, Exf, Eys,Eyf
	RANGE DEDA
	NONSPECIFIC_CURRENT i_m
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(um) = (micron)
} 

PARAMETER {
	x = 0 (um) : spatial coords  in microns
	y = 0 (um)
	z = 0 (um)
	XX = 0 (um)
	YY = 0 (um)
	FX = 0 (um)
	FY = 0 (um)
	DX=0 (um)
	DY=0 (um)
	Exs=0  (V millisec/m amp)
	Exf=0  (V millisec/m amp)
	Eys=0  (V millisec/m amp)
	Eyf=0  (V millisec/m amp)
	DEDA=0  (V millisec/m um amp)
}

ASSIGNED {
	rx  (millisec/m2)
	is  (amp/millisec)
	i_m (milliamp/cm2)
}

INITIAL {
	i_m=0
}

BREAKPOINT{
	i_m=-(0.1)*rx*is
}
