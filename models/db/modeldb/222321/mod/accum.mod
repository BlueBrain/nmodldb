TITLE accum :Ion accumulation, potassium buffering and diffusion between compartments and to the bath

NEURON {
	SUFFIX accum
	USEION na READ ina, nao, nai WRITE nai, nao, ina
	USEION k  READ ik, ko, ki WRITE ki, ko, ik
	USEION ca READ ica, cai WRITE cai, cao, ica VALENCE 2
	USEION cl READ icl, cli, clo WRITE cli, clo, icl VALENCE -1
	
	RANGE volin, volout, volpyram
	RANGE nap, kp, cap, clp
	POINTER diamp, inap, ikp, icap, iclp
	GLOBAL setnap, setkp, setcap, setclp, TotalBuffer
}

UNITS	{
	(um3)		= (liter/1e15)
	(mV)		= (millivolt)
	(mA)		= (milliamp)
	FARADAY		= 96485.309 (coul/mole)
	(molar)		= (1/liter)
	(mM)		= (millimolar)
	PI		= (pi) (1)
	R 		= (k-mole) (joule/degC)
}

PARAMETER {
	Difna = 0.1  (um2/ms)
	Difk  = 1.96	(um2/ms)
	Difca = 0.6	(um2/ms)
	Difcl = 2.03    (um2/ms)
	
	k1buf = 0.0008 (/ms)
	TotalBuffer = 500 (mM)
	x0 = 15 (mM)          

	nab = 140 (mM) 
	kb  = 3.5 (mM)
	clb = 135 (mM)
	
	s = 2e4

	setvolin=1
	setvolout=1
	setvolpyram=1

	setnap
	setkp
	setcap
	setclp
}

ASSIGNED {
	ina		(mA/cm2)
	ik		(mA/cm2)
	ica		(mA/cm2)
	icl		(mA/cm2)
	
	inap		(mA/cm2)
	ikp		(mA/cm2)
	icap		(mA/cm2)
	iclp		(mA/cm2)
	
	diam		(um)
	diamp		(um)
	
	naflux[3]
	kflux[3]
	caflux[3]
	clflux[3]
	
	dif[3]
	
	nai ki cai cli	:Concentrations in interneuron
	nao ko cao clo	:Concentrations in extracellular space
  	nap kp cap clp	:Concentrations in pyramidal neuron

	volin     	:Interneuron volume
	volout		:Extracellular volume
	volpyram	:Pyramidal neuron volume

	Kd (/mM)
	B0 (mM)
}

STATE { na[3] k[3] ca[3] cl[3] (mM) <1e-4> kbuf Buffer KBuffer (mM) vol[3] }

BREAKPOINT {
	SOLVE state METHOD sparse
}

INITIAL {
	na[0]=nai
	na[1]=nao
	na[2]=setnap
	k[0]=ki
	k[1]=ko
	k[2]=setkp
	ca[0]=cai
	ca[1]=cao
	ca[2]=setcap
	cl[0]=cli
	cl[1]=clo
	cl[2]=setclp
	
	:Potassium buffer
	Kd = k2buf(k[1])/k1buf
        kbuf = 0
	B0 = TotalBuffer/(1+Kd*k[1])
	Buffer = B0
	KBuffer = TotalBuffer - B0
	
	vol[0]=setvolin
	vol[1]=setvolout
	vol[2]=setvolpyram

	volin=vol[0]
	volout=vol[1]
	volpyram=vol[2]

	nap=na[2]
	kp=k[2]
	cap=ca[2]
	clp=cl[2]
}

KINETIC state {

	COMPARTMENT i, vol[i]*diam*diam*PI/4 { na k ca cl }

	LONGITUDINAL_DIFFUSION i, Difna*diam*diam*PI*vol[i]/4 { na }
	LONGITUDINAL_DIFFUSION i, Difk*diam*diam*PI*vol[i]/4  { k  }
	LONGITUDINAL_DIFFUSION i, Difca*diam*diam*PI*vol[i]/4 { ca }
	LONGITUDINAL_DIFFUSION i, Difcl*diam*diam*PI*vol[i]/4 { cl }

	:INTRACELLULAR INTERNEURON
	naflux[0] = -(ina*diam)*PI*(1e4)/FARADAY
        kflux[0]  = -(ik*diam)*PI*(1e4)/FARADAY
	caflux[0] = -(ica*diam)*PI*(1e4)/(FARADAY*2)
	clflux[0] = -(icl*diam)*PI*(1e4)/(FARADAY*-1)
	
	:COMMON EXTRACELLULAR SPACE	
	naflux[1] = (ina*diam+inap*diamp)*PI*(1e4)/FARADAY
        kflux[1]  = (ik*diam+ikp*diamp)*PI*(1e4)/FARADAY
	caflux[1] = (ica*diam+icap*diamp)*PI*(1e4)/(FARADAY*2)
	clflux[1] = (icl*diam+iclp*diamp)*PI*(1e4)/(FARADAY*-1)
	
	:INTRACELLULAR PYRAMIDAL
	naflux[2] = -(inap*diamp)*PI*(1e4)/FARADAY
	kflux[2]  = -(ikp*diamp)*PI*(1e4)/FARADAY
	caflux[2] = -(icap*diamp)*PI*(1e4)/(FARADAY*2)
	clflux[2] = -(iclp*diamp)*PI*(1e4)/(FARADAY*-1)
	
	:RADIAL DIFFUSION (TO THE BATH)
	dif[1]    = (Difna*geom(diam)*(nab-na[1])/s)*(setvolout*diam*diam*PI/4)
	dif[2]    = (Difk*geom(diam)*(kb-k[1]-kbuf)/s)*(setvolout*diam*diam*PI/4)
	dif[3]    = (Difcl*geom(diam)*(clb-cl[1])/s)*(setvolout*diam*diam*PI/4)

	~  na[0] << (naflux[0])
        ~  k[0]  << (kflux[0])
	~  ca[0] << (caflux[0])
	~  cl[0] << (clflux[0])
	~  na[1] << (naflux[1]+dif[1])
	~  k[1]  << (kflux[1]+dif[2])
	~  kbuf  << (-k2buf(k[1]+kbuf)*(k[1]+kbuf)*Buffer+k1buf*KBuffer)
  	~  Buffer << (-k2buf(k[1]+kbuf)*(k[1]+kbuf)*Buffer+k1buf*KBuffer)
  	~  KBuffer << (k2buf(k[1]+kbuf)*(k[1]+kbuf)*Buffer-k1buf*KBuffer)
	~  ca[1] << (caflux[1])
	~  cl[1] << (clflux[1]+dif[3])
	~  na[2] << (naflux[2])
	~  k[2]  << (kflux[2])
	~  ca[2] << (caflux[2])
	~  cl[2] << (clflux[2])

	volin=vol[0]
	volout=vol[1]
	volpyram=vol[2]

	nai=na[0] nao=na[1] nap=na[2]
	 ki=k[0]   ko=k[1]+kbuf   kp=k[2]
	cai=ca[0] cao=ca[1] cap=ca[2]
	cli=cl[0] clo=cl[1] clp=cl[2]
}

FUNCTION k2buf(x(mM)) {
	TABLE FROM -150 TO 150 WITH 200
	k2buf = k1buf/(1+exp((x-x0)/(-1.09)))
}

FUNCTION dr(x(um)) {
	TABLE FROM 0 TO 15 WITH 10
	dr = x*(sqrt(setvolout+1)-1)/2
}

FUNCTION geom(x(um)) {
	TABLE FROM 0 TO 15 WITH 10
	geom = PI*(x+2*dr(x))/(dr(x/2)*setvolout*x*x*PI/4)
}
