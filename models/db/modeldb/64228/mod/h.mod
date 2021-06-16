: Rod PhotoReceptor h channel

NEURON 
{
	SUFFIX h
	
	
	NONSPECIFIC_CURRENT ih
	RANGE  gh, ghbar,  Vhalfh, Sh
	RANGE  eh, aoh
	
	

}

UNITS
{
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millimho)
}

PARAMETER
{
	
	ghbar = 1.41 (mS/cm2) <0,1e9>
	
        eh  = -32 (mV)
       
        Vhalfh=-82  (mV)
       
        Sh = -5.33    (mV)
       
	aoh  = 1.64 (/s)
	
	
	

}

STATE
{
	
	nh
	
}

ASSIGNED
{
	v (mV)
	
	ih (mA/cm2)
	
	infh
	tauh   (ms)
	
	gh (mho/cm2)
	


}

INITIAL
{
	rate(v)
	
	nh  = infh
	
}

BREAKPOINT
{
	SOLVE states METHOD cnexp
	
	gh  = (0.001)*ghbar*nh
	ih  = gh*(v - eh)
	
	
}

DERIVATIVE states
{
	rate(v)
	
	nh'  = (infh  - nh )/tauh
	
}




FUNCTION alphah(v(mV)) (/ms)
{ 
	alphah = 0.001*aoh*exp( (v - Vhalfh)/(2*Sh)   )

}

FUNCTION betah(v(mV)) (/ms)
{ 
	betah =0.001*aoh*exp( - ( v-Vhalfh)/(2*Sh) )

}




PROCEDURE rate(v (mV))
{
        LOCAL  ah, bh

	
	
	ah = alphah(v)
	bh = betah(v)
	tauh = 1/(ah + bh)
	infh = ah*tauh
	
	

}

