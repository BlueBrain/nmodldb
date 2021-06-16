    COMMENT
    Point process for generation of multiple step currents at different amplitudes over time.
    ENDCOMMENT

    NEURON {
        POINT_PROCESS mStepIClamp
        RANGE del, dur
		RANGE ip1, ip2, ip3, ip4, ip5
		RANGE pon1, pon2, pon3, pon4, pon5
		RANGE poff1, poff2, poff3, poff4, poff5
		RANGE i0, s
        ELECTRODE_CURRENT i
    }

    UNITS {
        (nA) = (nanoamp)
    }

    PARAMETER {
        del=0   (ms)
        dur=10000   (ms)
    		
		ip1=0.5 (nA)
		ip2=-0.3 (nA)
		ip3=0.27 (nA)
		ip4=0.5 (nA)
		ip5=-0.37 (nA)
			
		pon1=1000 (ms)
		pon2=3000 (ms)
		pon3=5000 (ms)
		pon4=6500 (ms)
		pon5=9000 (ms)
			
		poff1=1300 (ms)
		poff2=3300 (ms)
		poff3=9000 (ms)
		poff4=6800 (ms)
		poff5=9300 (ms)
			
		i0=6.78 (nA)
		s=17			
    }

    ASSIGNED {
        i (nA)
    }

    BREAKPOINT {
        at_time(del)
        at_time(del + dur)
		
		i = i0 + s*((heav(poff1-t)*heav(t-pon1)*ip1) + (heav(poff2-t)*heav(t-pon2)*ip2-0.01) + (heav(poff3-t)*heav(t-pon3)*ip3) + (heav(poff4-t)*heav(t-pon4)*ip4) + (heav(poff5-t)*heav(t-pon5)*ip5))
    }
	
	FUNCTION heav (x) {
		if (x < 0) {heav = 0}
		else {heav = 1}
	}