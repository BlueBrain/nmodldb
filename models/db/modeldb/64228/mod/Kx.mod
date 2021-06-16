: Rod PhotoReceptor Kx channel

NEURON 
{
	SUFFIX Kx
	
	USEION Kx WRITE iKx VALENCE 1
	
	
	RANGE gKx, gKxbar, VhalfKx, SKx
	RANGE eKx, aoKx
	
	

}

UNITS
{
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millimho)
}

PARAMETER
{
	gKxbar = 0.85 (mS/cm2) <0,1e9>
	
        eKx = -74 (mV)
                
        VhalfKx=-49.9 (mV)
        
        
	aoKx = 0.66 (/s)
	
	SKx = 5.7     (mV) 
	
   	

}

STATE
{
	nKx
	
}

ASSIGNED
{
	v (mV)
	
	
	iKx (mA/cm2)
	
	infKx
	tauKx  (ms)
	
	gKx (mho/cm2)
	
}

INITIAL
{
	rate(v)
	nKx = infKx
	
}

BREAKPOINT
{
	SOLVE states METHOD cnexp
	gKx = (0.001)*gKxbar*nKx
	iKx = gKx*(v - eKx)
	
	
}

DERIVATIVE states
{
	rate(v)
	nKx' = (infKx - nKx)/tauKx
	
}


FUNCTION alphaKx(v(mV)) (/ms)
{ 
	alphaKx = 0.001*aoKx * exp( (v - VhalfKx)/(2*SKx) )
}

FUNCTION betaKx(v(mV))  (/ms)
{ 
	betaKx = 0.001*aoKx*exp( - ( v-VhalfKx)/(2*SKx) )
}


PROCEDURE rate(v (mV))
{
        LOCAL aKx, bKx

	aKx = alphaKx(v)
	bKx = betaKx(v)
	tauKx = 1/(aKx + bKx)
	infKx = aKx*tauKx
	
}

