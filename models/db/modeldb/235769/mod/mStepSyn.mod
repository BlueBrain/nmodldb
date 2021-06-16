    COMMENT
    Point process for generation of multiple step changes in synaptic conductance over time.
	ENDCOMMENT

    NEURON {
        POINT_PROCESS StepSyn
        RANGE del, dur
		RANGE ip1, ip2, ip3, ip4, ip5
		RANGE pon1, pon2, pon3, pon4, pon5
		RANGE poff1, poff2, poff3, poff4, poff5
		RANGE i0, s, g
        :ELECTRODE_CURRENT i
		NONSPECIFIC_CURRENT i
    }

    UNITS {
        (nA) = (nanoamp)
		(mV) = (millivolt)
		(uS) = (microsiemens)
    }

    PARAMETER {
		e=0	(mV)
        del=0   (ms)
        dur=10000   (ms)
    		
		ip1=0.6 (uS)	<0,1e9>
		ip2=-0.32 (uS)	<0,1e9>
		ip3=0.07 (uS)	<0,1e9>
		ip4=0.15 (uS)	<0,1e9>
		ip5=-0.37 (uS)	<0,1e9>
			
		pon1=1000 (ms)
		pon2=3000 (ms)
		pon3=5000 (ms)
		pon4=7000 (ms)
		pon5=9000 (ms)
			
		poff1=1200 (ms)
		poff2=3200 (ms)
		poff3=9000 (ms)
		poff4=7200 (ms)
		poff5=9200 (ms)
			
		i0=6.5 (uS)	<0,1e9>
		s=17			
    }

    ASSIGNED {
        v (mV) 
		i (nA)  
		g (uS)
    }

    BREAKPOINT {
        at_time(del)
        at_time(del + dur)
		
		g = i0 + s*((heav(poff1-t)*heav(t-pon1)*ip1) + (heav(poff2-t)*heav(t-pon2)*ip2) + (heav(poff3-t)*heav(t-pon3)*ip3) + (heav(poff4-t)*heav(t-pon4)*ip4) + (heav(poff5-t)*heav(t-pon5)*ip5))
		i = g*(v - e)
	}
	
	FUNCTION heav (x) {
		if (x < 0) {heav = 0}
		else {heav = 1}
	}