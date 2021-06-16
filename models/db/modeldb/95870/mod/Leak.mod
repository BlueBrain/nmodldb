: Rod  Photoreceptor Leakage

NEURON 
{
	SUFFIX Leak
	

	
	NONSPECIFIC_CURRENT il
	
	
          
	RANGE glbar, el
	
	

	

}

UNITS
{
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millimho)
        
}	

PARAMETER
{
       
       
       : leak
        glbar=0.6 (mS/cm2)
       
	: for 100uM^2 area, =1nS/1GOhm
	 el=-55 (mV)

}


ASSIGNED
{
	v (mV)
	gl (mho/cm2)
	il  (mA/cm2)
          
}





BREAKPOINT
{
        gl= 1e-3*glbar :gl is 3 orders different from glbar

	il  = gl*(v-el)

	: the current is in the unit of mA/cm2
	
}

