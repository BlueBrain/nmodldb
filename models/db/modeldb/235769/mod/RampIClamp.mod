    COMMENT
    Point process for generation of ascending and descending current injection over time.
    ENDCOMMENT

    NEURON {
            POINT_PROCESS RampIClamp
            RANGE del, dur, pkamp, bias
            ELECTRODE_CURRENT i
    }

    UNITS {
            (nA) = (nanoamp)
          }

    PARAMETER {
            del=0   (ms)
            dur=10000   (ms)
            pkamp=20 (nA)
			bias=0 (nA)
    }

    ASSIGNED {
            i (nA)
    }

    BREAKPOINT {
           at_time(del)
           at_time(del + dur)

           if (t < del) {
          i=0   
       }else{ 
                if (t < del+dur/2) {
               i = (2*pkamp/dur)*t + bias
          }else{ 
               i = -(2*pkamp/dur)*(t-dur) + bias
    }}}