TITLE kca.mod   Calcium activated K channel

COMMENT

    Calcium activated Potassium channel

    Simplifictaion: Spike dependent current only!

    Created by Christian Roessert

ENDCOMMENT




UNITS {
        (nA) = (milliamp)
        (mV) = (millivolt)
}

NEURON {
  POINT_PROCESS KCa
  NONSPECIFIC_CURRENT ik
  RANGE dgkbar, egk, ctau, thresh, gk
}


PARAMETER {
  	dgkbar = 0.1 (mS/cm2) <0,1e9>
	egk   = -90 (mV)   
  	ctau =  50   (ms)
  	thresh = -20 (mV)
}


STATE {
        c
}

ASSIGNED {
  	v (mV)
	area (um2) : area of current segment (automatically available within NMODL, like v)
 	dgk (uS) : set to dgkbar * segment area at initialization
  	gk (uS) : will be 0 or dgk*c depending on recent spiking history
	ik (nA)
}

BREAKPOINT {
	SOLVE state METHOD cnexp 
	gk = dgk*c
	ik = gk*(v - egk) 

}

INITIAL {
  	dgk = dgkbar*area*1e-5 : because area will be in um2, but dgk is in uS and dgkbar in mS/cm2
   	gk = 0 : because at t = 0 we assume that the cell has not yet spiked
	c = 0
	net_send(0, 1)
}

DERIVATIVE state {     : exact when v held constant; integrates over dt step
	c' = -c/ctau
}

NET_RECEIVE (null) {
	if (flag==1) { : no spike has occured
		WATCH (v > thresh) 2 : detect spike
	}

	if (flag==2) { : spike has occured
		WATCH (v < thresh) 3 : detect spike off
        :c = c+1 : increase c and wait for next spike
		:net_send(0, 1)
   	}
	
	if (flag==3) { : 
		c = c+1 : increase c and wait for next spike
        at_time(t)
		net_send(0, 1)
   	}
}



