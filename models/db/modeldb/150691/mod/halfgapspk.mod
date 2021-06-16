NEURON {
	POINT_PROCESS HalfGapSpk
	: gap changes its conductance when it receives a NetCon event
	: it is critical that both sides of the gaps receive the same
	: events at the same time.
	ELECTRODE_CURRENT i
	RANGE g, i, vgap, meang, meant, rg, rt, drift
	THREADSAFE : Only true if every instance has its own distinct Random
	POINTER donotuse : A Normal Random generator with mean 1 and var 1.
	RANGE id : For polarity of rectification and testing.
		 : Should be equal and opposite for corresponding HalfGap
		 : and otherwise distinct. For proper simulation results,
		 : corresponding gaps should always have the same value
		 : of g.
	RANGE gmax, gmin, vhalf : Sigmoidal voltage sensitive conductance
		: parameters. See gv(x) below. The sign of id defines
		: the voltage polarity. If gmax == gmin, the gap is linear
		: and id is not used.
		: in pargap, gmin==gmax (linear) unless gmin is 0.
}

PARAMETER {
	gmax = 1 (nanosiemens)
	gmin = 1 (nanosiemens)
	vhalf = 0 (millivolt)
	slope4 = 10 (/millivolt)
	meang = 30 (nanosiemens)
	meant = 1000000 (ms)
	drift = 0
	rg=0
	rt=0
	id = 0
}

ASSIGNED {
	g (nanosiemens)
	v (millivolt)
	vgap (millivolt)
	i (nanoamp)
	donotuse
}

INITIAL {
}


: voltage sensitve gap conductance
: for global variable time step, should be continuous to high order so
: that performance does not suffer.
: Argument is relative voltage at the positive polarity side.
FUNCTION gv(x(millivolt))(nanosiemens) {
	: sigmoid x >> vhalf means gv = gmax, x << vhalf means g = gmin
	gv = (gmax - gmin)/(1 + exp(slope4*(vhalf - x))) + gmin
}

BREAKPOINT {
	LOCAL x
	if (gmax == gmin) { :linear gap junction
		g = gmax
		i = g * (vgap - v) * (.001)
	}else{
		: vgap > v means current is outward from this gap
		if (id > 0 ) {
			x = v - vgap :voltage relative to - side of gap
		}else if (id < 0){
			x = vgap - v : voltage relative to - side of gap
		}else{
VERBATIM
			assert(0);
ENDVERBATIM

		}
		g = gv(x)
		i = g * (vgap - v) * (.001)
	}
}

PROCEDURE getpar() {
	gmax=mynormrand(meang/1(nanosiemens),rg)*1(nanosiemens)
	if (gmax<0) {gmax=0}
	if (gmin != 0) {
		gmin = gmax
	}
	meang=meang+drift*meang
	rg=rg+drift*rg
}

NET_RECEIVE (w) {
	getpar()		:sets gmax,=gmin and next change
}

:Separate independent but reproducible streams for each instance.
:For proper functioning, it is important that hoc Random distribution be
: Random.Random123(id1, id2) <one could use MCellRan4 instead>
: Random.normal(1,1)
: and that corresponding HalfGap have the same id1, id2
: A condition for correctness, that can be tested from hoc, is that
: g (and also Random.seq()) for corresponding HalfGap have the same value.
: If this is the case, then simulations with different numbers of processes
: and different distibutions of gids should give quantitatively identical
: results with the fixed step method and  (if cvode.use_long_double(1))
: with the global variable time step method.

VERBATIM
double nrn_random_pick(void* r); 
void* nrn_random_arg(int argpos);
ENDVERBATIM


FUNCTION mynormrand(mean, var) {
VERBATIM
	if (_p_donotuse) {
		double x = nrn_random_pick(_p_donotuse);
		_lmynormrand = x*_lvar + _lmean;
	}else{
		_lmynormrand = _lmean;
	}
ENDVERBATIM
}

PROCEDURE setRandom() {
VERBATIM
 {
        void** pv = (void**)(&_p_donotuse);
        if (ifarg(1)) {
                *pv = nrn_random_arg(1);
        }else{
                *pv = (void*)0;
        }
 }
ENDVERBATIM
}

