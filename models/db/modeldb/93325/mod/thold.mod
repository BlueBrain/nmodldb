: Threshold Function, Bob Calin-Jageman 9/20/2002

: This model implements a threshold function described by 
: Getting, 1989 and utilized by Lieb and Frost, 1998.
: 
: 
: Created by Bob Calin-Jageman
: 
: Created 	9/20/2002
: Modified 	9/20/2002
:
: Mathcheck  - Complete    9/31/2002
: Unitscheck - Complete	   9/31/2002 
: 
: Explanation
: With every action potential (v.pre > threshold), the threshold
: is reset to RESET and then decays back to STEADYSTATE with time constant 
: DECAYTC.
: 
: This model is used in conjunction with a netcon object, which detects when
: the membrane potential of a neuron exceed a set threshold.  THe pointer
: nc_thresh is set to point at the threshold of the netcon.  Thus, this model
: causes a typical netcon to have a variable threshold that changes with
: activity.
: 
: This model must also receive input from the netcon to detect the occurance
: of spike.  The hoc code for this model should look like this
: 
:Example Usage
:
:		soma {sthold = new thold(0.5)}
:		sthold.steadystate = sthold_steadystate
:		sthold.reset = sthold_reset
:		sthold.decaytc = sthold_decaytc
:
:		soma stholdnc = new NetCon(&v(0.5), sthold, 0, 0, 1)
:		setpointer sthold.nc_thresh, stholdnc.threshold
: 
: References
: 	Getting, P.A. (1989) "Reconstruction of small neural networks" in 
: Methods in Neuronal Modeling: From Synapses to Networks (1st ed), Kock & Segev
: eds, MIT Press.
:
: 	Lieb JR & Frost WN (1997) "Realistic Simulation of the Aplysia Siphon
: Withdrawal Reflex Circuit: Roles of Circuit Elements in Producing Motor Output"
: p. 1249 */
: 

NEURON {
	POINT_PROCESS thold
	RANGE steadystate, reset, decaytc, lastspike, thresh, spikecount, prior2spike, burstmaxint, burstminsize, bursting, burstcount, burststart, burstno, freq
	POINTER nc_thresh
	}

UNITS {
	(mV) = (millivolt)
}

PARAMETER  {
	steadystate = -37.8 (mV)  : resting threshold level 
	reset = 200 (mV)	  : with each spike, threshold reset to here 
	decaytc = 11.5 (ms)	  : threshold decays back to steadystate with this time constant 
	nc_thresh (mV)		  : pointer.  Modifies the threshold for a netcon object 
	lastspike (ms)		  : time of last spike.  initialized to 0 don't modify extern
	thresh (mV)		  : current threshold 
	spikecount		  :
	prior2spike		  : time of spike before lasrt.  initialized to -1 don't modify extern
	burstmaxint = 1000
	burstminsize = 3
	bursting = 0
	burstcount = 0
	burststart = 0
	burstno = 0
	freq = 0
	}

ASSIGNED {
	v	(mV)
	}


BREAKPOINT {	
	thresh = steadystate + (reset - steadystate) * (exp((lastspike - t) / decaytc))
	nc_thresh = thresh  : set netcon object to reflect current thresh
}	

INITIAL {
	thresh = reset 
	nc_thresh = thresh
	lastspike = 0
	spikecount = 0
	prior2spike = -1
	bursting = 0
	burstcount = 0
	burststart = 0
	burstno = 0
	freq = 0
}


NET_RECEIVE(weight (microsiemens)) {
		if (spikecount > 0) {
		freq = (1/ (t - lastspike)* 1000) 
		if ((t - lastspike) < burstmaxint) {
			burstcount = burstcount + 1
			if (burstcount == 1) { burststart = lastspike }
			if (burstcount > burstminsize) { bursting = 1 }
		} 
		
		if ((t - lastspike) > burstmaxint) {
			bursting = 0
			burstcount = 0
			burststart = 0
		}
		}

		prior2spike = lastspike
		lastspike = t     : recieved from netcon when v > thresh.  
		spikecount = spikecount + 1  : increment spike count
		
	
}