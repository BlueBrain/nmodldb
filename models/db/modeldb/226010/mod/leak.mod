TITLE  leak current

COMMENT
written for NEURON by Antonios Dougalis, 25 Feb 2015, London
based on voltage clamp data from Dougalis et al., 2017 J Compu Neurosci 
ENDCOMMENT

UNITS {
        (S) = (siemens)
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX leak
        NONSPECIFIC_CURRENT ileak 
        RANGE gleakbar,ileak,eleak

}
 
PARAMETER {
        v   (mV)
        dt  (ms)
		gleakbar  = 0.00002  (S/cm2)
        eleak = -55 (mV)
		
}
 
 
ASSIGNED {
        ileak (mA/cm2)
        
}
 
BREAKPOINT {
         ileak = gleakbar*(v - eleak)  

		
}
 
 
