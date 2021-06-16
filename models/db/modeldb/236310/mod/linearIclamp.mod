COMMENT
Since this is an electrode current, positive values of i depolarize the cell
and in the presence of the extracellular mechanism there will be a change
in vext since i is not a transmembrane current but a current injected
directly to the inside of the cell.
ENDCOMMENT

NEURON {
        POINT_PROCESS linearIclamp
        RANGE del, dur, pkamp, freq, phase, bias
        ELECTRODE_CURRENT i
}

UNITS {
        (nA) = (nanoamp)
             }

PARAMETER {
        del=215   (ms)
        dur=1000   (ms)
        slamp=0.3 (nA)
}

ASSIGNED {
        i (nA)
}

BREAKPOINT {
       if (t < del) {
      i=0   
   }else{  
            if (t < del+dur) {
           i = slamp*(t-del)
      }else{  
           i = 0
}}}
