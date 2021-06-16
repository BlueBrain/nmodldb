:Comment :
:Reference : :		Reuveni, Friedman, Amitai, and Gutnick, J.Neurosci. 1993

NEURON	{
	SUFFIX Ca_HVA
	USEION ca READ eca,cai WRITE ica
	RANGE gCa_HVAbar, gCa, ica 
	RANGE ma,mb,ha,hb
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gCa_HVAbar = 0.00001 (S/cm2) 
	ma = 1
	mb = 1
	ha = 1
	hb = 1
	
}

ASSIGNED	{
	v	(mV)
	eca	(mV)
	ica	(mA/cm2)
	gCa	(S/cm2)
	cai (mM)

	mInf
	mTau
	mAlpha
	mBeta
	hInf
	hTau
	hAlpha
	hBeta
	zInf
}

STATE	{ 
	m
	h
	z
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gCa = gCa_HVAbar*m*m*h
	ica = gCa*(v-eca)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
}

INITIAL{
	rates()
	m = mInf
	h = hInf
}

PROCEDURE rates(){
	UNITSOFF
        if((v == -27 ||v == -75 ||v == -13  ||v == -15  ) ){        
            v = v + 0.0001 
        }
		mAlpha =  ma*(0.055*(-27-v))/(exp((-27-v)/3.8) - 1)        
		mBeta  =  mb*(0.94*exp((-75-v)/17))
		mInf = mAlpha/(mAlpha + mBeta)
		mTau = 1/(mAlpha + mBeta)

		hAlpha =  ha*(0.000457*exp((-13 - v )/50))
		hBeta  =  hb*(0.0065/(exp((-v - 15)/28)+1))
		hInf = hAlpha/(hAlpha + hBeta)
		hTau = 1/(hAlpha + hBeta)
	UNITSON
}
