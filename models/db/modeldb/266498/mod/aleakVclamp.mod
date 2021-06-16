: A passive leak current
NEURON {
	SUFFIX aleak
	USEION na READ ena WRITE ina
	USEION ca READ cai, cao WRITE ica
	RANGE i, ina, ica, gbna, gbca
}
UNITS {
	(S)=(siemens)
	(mV)=(millivolt)
	(mA)=(milliamp)
	F = (faraday) (coulombs)
	(molar) = (1/liter)
	(mM)    = (millimolar)
}
PARAMETER {
	gbna=1.14945E-05 (S/cm2) <0, 1e9>
	gbca=3.00626E-06 (S/cm2) <0, 1e9>
	
	R=8.314 (joule/degC): Gas Constant
	z=2 : Charge of Ca ion
	ecaoffset=78.7 (mV)
}
ASSIGNED {
	v (mV)
	ena (mV)
	eca (mV)
	ica (mA/cm2)
	ina (mA/cm2)
	celsius (degC)
	ecaleak	(mV)
	cai(mM)
	cao (mM)
	
}
BREAKPOINT { 

	ecaleak=(1000)*(R*(celsius+273.15)/z/F*log(cao/cai))-ecaoffset : Equation for eca given in Schild 1994
	ina = gbna*(v - ena) 
	ica = gbca*(v - ecaleak) 

}