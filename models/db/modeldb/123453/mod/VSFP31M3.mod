TITLE Voltage sensor protein VSFP3.1

COMMENT
Model calculates the displacement current and fluorescence response produced by the voltage-sensing domain in VSFP3.1
Kinetic parameters from least-square fits to experimental data provided by Alica Lundby of VSFP3.1 expressed in PC12 cells

Version 2.0

KINETIC SCHEME: 
8-state-Markov-Process

Sensor states: closed = S1n, S2n; open = S1p, S2p
Reporter: off = Rn; on = Rp

Reference: Akemann et al. Biophys. J. (2009) 96: 3959-3976

Laboratory for Neuronal Circuit Dynamics
RIKEN Brain Science Institute, Wako City, Japan
http://www.neurodynamics.brain.riken.jp

Date of Implementation: December 2008
Contact: akemann@brain.riken.jp

ENDCOMMENT

NEURON {
	SUFFIX VSFP31M3
	NONSPECIFIC_CURRENT i
	GLOBAL S1ONzero, S1OFFzero, S2niONzero, S2niOFFzero, S2ipONzero, S2ipONzero, S12pONzero, S12pOFFzero, S12nONzero, S12nOFFzero
	GLOBAL R1nONzero, R1nOFFzero, R1pONzero, R1pOFFzero, R2pONzero, R2pOFFzero, R2nONzero, R2nOFFzero
	GLOBAL zGateS1, zGateS2, zGateS12p, zGateS12n
	GLOBAL deltaGateS1, deltaGateS2ni, deltaGateS2i, deltaGateS2ip, deltaGateS12p, deltaGateS12n
	GLOBAL deltaF, Fhalf
	GLOBAL baseline
	RANGE nc
	RANGE fluoSignal
	RANGE fluoActivation
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(nA) = (nanoamp)
	(pA) = (picoamp)
	(S)  = (siemens)
	(mS) = (millisiemens)
	(nS) = (nanosiemens)
	(pS) = (picosiemens)
	(um) = (micron)
	(molar) = (1/liter)
	(mM) = (millimolar)		
}

CONSTANT {
	e0 = 1.60217646e-19 (coulombs)
	kB = 1.3806505e-23 (joule/kelvin) 
	q10Gate = 1.43 (1)			: q10 = 0.9 (70 mV); 0.6 (-30 mV)  ON
	q10Fluo = 1.67 (1) 				
	tempGate = 25 (degC)
	tempFluo = 25 (degC)
}

PARAMETER {
	baseline = 1 (1)	: 0 = fluorescence baseline set to 0; 

	nc = 0 (1/cm2)
	
	S1ONzero = 0.48 (1/ms)
	S1OFFzero = 0.074 (1/ms)

	S2niONzero = 0.4 (1/ms)
	S2niOFFzero = 0.38 (1/ms)

	S2ipONzero = 2 (1/ms)
	S2ipOFFzero = 0.1 (1/ms)

	S12pONzero = 0.014 (1/ms)
	S12pOFFzero = 0.0066 (1/ms)

	S12nONzero = 1e-9 (1/ms)
	S12nOFFzero = 0.002 (1/ms)

	zGateS1 = 1.2 (1)
	zGateS2 = 1.2 (1)
	zGateS12p = 0.3 (1)
	zGateS12n = 0.3 (1)
	
	deltaGateS1 = 0.35 (1)		: location of the S1 transition state (between 0 = internal side to 1 = external side)
	deltaGateS2i = 0.2 (1)		: location of the transition state of the intermediate state (Si) in S2
	deltaGateS2ni = 0.15 (1)	: location of the transition state in the reaction of S2n to S2i
	deltaGateS2ip = 0.3 (1)		: location of the transition state in the reaction of S2i to S2p
	deltaGateS12p = 0.4 (1)		: location of the transition state in the reaction of S1p to S2p
	deltaGateS12n = 0.2 (1)		: location of the transition state in the reaction of S1n to S2n

	R1nONzero = 1e-12 (1/ms)
	R1nOFFzero = 2 (1/ms)
	
	R1pONzero = 1 (1/ms)
	R1pOFFzero = 0.7 (1/ms)

	R2pONzero = 2 (1/ms)
	R2pOFFzero = 1e-12 (1/ms)
	
	R2nONzero = 1e-9 (1/ms)
	R2nOFFzero = 0.028 (1/ms)

	deltaF = -0.0105 (1)			: maximum fluoresence modulation
	Fhalf = 1 (1)				: fluorescence ratio at vhalf
}

ASSIGNED {
	celsius (degC)
	v (mV)
	
	i (mA/cm2)
	
	S1ON (1/ms)
	S1OFF (1/ms)

	S2niON (1/ms)
	S2niOFF (1/ms)
	S2ipON (1/ms)
	S2ipOFF (1/ms)

	S12pON (1/ms)
	S12pOFF (1/ms)

	S12nON (1/ms)
	S12nOFF (1/ms)

	R1nON (1/ms)
	R1nOFF (1/ms)
	
	R1pON (1/ms)
	R1pOFF (1/ms)

	R2pON (1/ms)
	R2pOFF (1/ms)

	R2nON (1/ms)
	R2nOFF (1/ms)

	qtGate (1)
	qtFluo (1)

	fluoSignal (1)				: Fluorescence response
	fluoActivation (1)	
	fluoInit (1)
}

STATE {
	S1nRn FROM 0 TO 1
	S1pRn FROM 0 TO 1
	S2pRn FROM 0 TO 1
	S2nRn FROM 0 TO 1
	S2iRn FROM 0 TO 1	

	S1nRp FROM 0 TO 1
	S1pRp FROM 0 TO 1
	S2pRp FROM 0 TO 1
	S2nRp FROM 0 TO 1
	S2iRp FROM 0 TO 1
}

INITIAL {
	qtGate = q10Gate^((celsius-tempGate)/10 (degC))
	qtFluo = q10Fluo^((celsius-tempFluo)/10 (degC))

	if ( deltaGateS2ni > deltaGateS2i ) { deltaGateS2ni = deltaGateS2i }
	if ( deltaGateS2ip < deltaGateS2i ) { deltaGateS2ip = deltaGateS2i }

	rateConst(v)
	SOLVE reaction STEADYSTATE sparse

	if ( baseline == 0 ) {
		fluoInit = Fhalf + deltaF * ( S1nRp + S1pRp + S2pRp + S2nRp + S2iRp - 0.5 )
		} else {
		fluoInit = 0
		}
}

BREAKPOINT {
	SOLVE reaction METHOD sparse
	i = nc * (1e6) * e0 * ( zGateS1 * gate1Flip() + zGateS2 * gate2Flip() + zGateS12p * gate12pFlip() + zGateS12n * gate12nFlip() )
 	
	fluoActivation = S1nRp + S1pRp + S2pRp + S2nRp + S2iRp 
	fluoSignal = Fhalf + deltaF * ( S1nRp + S1pRp + S2pRp + S2nRp + S2iRp - 0.5 ) - fluoInit
}

KINETIC reaction {
	rateConst(v)
	~ S1nRn <-> S1pRn		(S1ON, S1OFF)
	~ S2nRn <-> S2iRn		(S2niON, S2niOFF)
	~ S2iRn <-> S2pRn		(S2ipON, S2ipOFF)

	~ S1pRn <-> S2pRn		(S12pON, S12pOFF)
	~ S1nRn <-> S2nRn		(S12nON, S12nOFF)

	~ S1nRp <-> S1pRp		(S1ON, S1OFF)
	~ S2nRp <-> S2iRp		(S2niON, S2niOFF)
	~ S2iRp <-> S2pRp		(S2ipON, S2ipOFF)

	~ S1pRp <-> S2pRp		(S12pON, S12pOFF)
	~ S1nRp <-> S2nRp		(S12nON, S12nOFF)

	~ S1nRn <-> S1nRp		(R1nON, R1nOFF)
 	~ S1pRn <-> S1pRp		(R1pON, R1pOFF)
	~ S2pRn <-> S2pRp		(R2pON, R2pOFF)
	~ S2nRn <-> S2nRp		(R2nON, R2nOFF)


CONSERVE S1nRn + S1pRn + S2pRn + S2nRn + S2iRn + S1nRp + S1pRp + S2pRp + S2nRp + S2iRp = 1
}

PROCEDURE rateConst( v(mV) ) {
	
	S1ON = qtGate * S1ONzero * exp( zGateS1 * e0 * deltaGateS1 * v / ( (1000) * kB * celsiusTOkelvin( celsius ) ) )
	S1OFF = qtGate * S1OFFzero * exp( -zGateS1 * e0 * (1-deltaGateS1) * v / ( (1000) * kB * celsiusTOkelvin( celsius ) ) )
	
	S2niON = qtGate * S2niONzero * exp( zGateS2 * e0 * deltaGateS2ni * v / ( (1000) * kB * celsiusTOkelvin( celsius ) ) )
	S2niOFF = qtGate * S2niOFFzero * exp( -zGateS2 * e0 * (deltaGateS2i-deltaGateS2ni) * v / ( (1000) * kB * celsiusTOkelvin( celsius ) ) )

	S2ipON = qtGate * S2ipONzero * exp( zGateS2 * e0 * (deltaGateS2ip-deltaGateS2i) * v / ( (1000) * kB * celsiusTOkelvin( celsius ) ) )
	S2ipOFF = qtGate * S2ipOFFzero * exp( -zGateS2 * e0 * (1-deltaGateS2ip) * v / ( (1000) * kB * celsiusTOkelvin( celsius ) ) )

	S12pON = qtGate * S12pONzero * exp( zGateS12p * e0 * deltaGateS12p * v / ( (1000) * kB * celsiusTOkelvin( celsius ) ) )
	S12pOFF = qtGate * S12pOFFzero * exp( -zGateS12p * e0 * (1-deltaGateS12p) * v / ( (1000) * kB * celsiusTOkelvin( celsius ) ) )

	S12nON = qtGate * S12nONzero * exp( zGateS12n * e0 * deltaGateS12n * v / ( (1000) * kB * celsiusTOkelvin( celsius ) ) )
	S12nOFF = qtGate * S12nOFFzero * exp( -zGateS12n * e0 * (1-deltaGateS12n) * v / ( (1000) * kB * celsiusTOkelvin( celsius ) ) )

	R1nON = qtFluo * R1nONzero
	R1nOFF = qtFluo * R1nOFFzero

	R1pON = qtFluo * R1pONzero
	R1pOFF = qtFluo * R1pOFFzero

	R2pON = qtFluo * R2pONzero
	R2pOFF = qtFluo * R2pOFFzero

	R2nON = qtFluo * R2nONzero
	R2nOFF = qtFluo * R2nOFFzero
}		

FUNCTION gate1Flip() (1/ms) {
	gate1Flip = S1ON * ( S1nRn + S1nRp ) - S1OFF * ( S1pRn + S1pRp )
}

FUNCTION gate2Flip() (1/ms) {
	gate2Flip = deltaGateS2i * ( S2niON * ( S2nRn + S2nRp ) - S2niOFF * ( S2iRn + S2iRp ) ) + (1-deltaGateS2i) * ( S2ipON * ( S2iRn + S2iRp ) - S2ipOFF * ( S2pRn +S2pRp ) )
}

FUNCTION gate12pFlip() (1/ms) {
	gate12pFlip = S12pON * ( S1pRn + S1pRp ) - S12pOFF * ( S2pRn + S2pRp )
}

FUNCTION gate12nFlip() (1/ms) {
	gate12nFlip = S12nON * ( S1nRn + S1nRp ) - S12nOFF * ( S2nRn + S2nRp )
}

FUNCTION celsiusTOkelvin ( c (degC) ) (kelvin) {
UNITSOFF
	celsiusTOkelvin = 273.15 + c
UNITSON
}

