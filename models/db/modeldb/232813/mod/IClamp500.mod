: IClamp500.mod codes 500 Hz pulse current injection during current-clamp recoding.
:
: Takaki Watanabe
: wtakaki@m.u-tokyo.ac.jp

NEURON {
POINT_PROCESS IClamp500
RANGE del, dur, invl, amp, i 
ELECTRODE_CURRENT i
}

UNITS {
(nA) = (nanoamp)
   }

PARAMETER {
  	del (ms)
   	dur (ms)	<0,1e9>
	invl (ms) <0,1e9>
    amp (nA)
   }
ASSIGNED { i (nA) }

INITIAL {
   	i = 0
    }
  
 BREAKPOINT {
 	at_time(del)
 	at_time(del+dur)
	at_time(del+dur+invl)
   	if (t < del + dur  && t >= del ) {
   	i = amp
	}else if (t < del + dur + 1*(dur + invl) && t >= del + 1*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 2*(dur + invl) && t >= del + 2*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 3*(dur + invl) && t >= del + 3*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 4*(dur + invl) && t >= del + 4*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 5*(dur + invl) && t >= del + 5*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 6*(dur + invl) && t >= del + 6*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 7*(dur + invl) && t >= del + 7*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 8*(dur + invl) && t >= del + 8*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 9*(dur + invl) && t >= del + 9*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 10*(dur + invl) && t >= del + 10*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 11*(dur + invl) && t >= del + 11*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 12*(dur + invl) && t >= del + 12*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 13*(dur + invl) && t >= del + 13*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 14*(dur + invl) && t >= del + 14*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 15*(dur + invl) && t >= del + 15*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 16*(dur + invl) && t >= del + 16*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 17*(dur + invl) && t >= del + 17*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 18*(dur + invl) && t >= del + 18*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 19*(dur + invl) && t >= del + 19*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 20*(dur + invl) && t >= del + 20*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 21*(dur + invl) && t >= del + 21*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 22*(dur + invl) && t >= del + 22*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 23*(dur + invl) && t >= del + 23*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 24*(dur + invl) && t >= del + 24*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 25*(dur + invl) && t >= del + 25*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 26*(dur + invl) && t >= del + 26*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 27*(dur + invl) && t >= del + 27*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 28*(dur + invl) && t >= del + 28*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 29*(dur + invl) && t >= del + 29*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 30*(dur + invl) && t >= del + 30*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 31*(dur + invl) && t >= del + 31*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 32*(dur + invl) && t >= del + 32*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 33*(dur + invl) && t >= del + 33*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 34*(dur + invl) && t >= del + 34*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 35*(dur + invl) && t >= del + 35*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 36*(dur + invl) && t >= del + 36*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 37*(dur + invl) && t >= del + 37*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 38*(dur + invl) && t >= del + 38*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 39*(dur + invl) && t >= del + 39*(dur + invl)) {
	i = amp  
    }else if (t < del + dur + 40*(dur + invl) && t >= del + 40*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 41*(dur + invl) && t >= del + 41*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 42*(dur + invl) && t >= del + 42*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 43*(dur + invl) && t >= del + 43*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 44*(dur + invl) && t >= del + 44*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 45*(dur + invl) && t >= del + 45*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 46*(dur + invl) && t >= del + 46*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 47*(dur + invl) && t >= del + 47*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 48*(dur + invl) && t >= del + 48*(dur + invl)) {
	i = amp
	}else if (t < del + dur + 49*(dur + invl) && t >= del + 49*(dur + invl)) {
	i = amp 	
	}else{
	i = 0
   }
   }
