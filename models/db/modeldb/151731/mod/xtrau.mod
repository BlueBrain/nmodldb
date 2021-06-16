: xtrau.mod, based on xtra.mod,v 1.3 2009/02/24 00:52:07
: the "u" stands for "uniform" as in "uniform field"
: because this modification was made to accommodate a
: uniform extracellular field,
: as would occur between two parallel plates

COMMENT
This mechanism is intended to be used in conjunction 
with the extracellular mechanism.  Pointers specified 
at the hoc level must be used to connect the 
extracellular mechanism's e_extracellular and i_membrane 
to this mechanism's ex and im, respectively.

xtrau does four useful things:

1. Serves as a target for Vector.play() to facilitate 
extracellular stimulation.  Assumes that one has initialized 
a Vector to hold the time sequence of the extracellular field
strength E.

This Vector is to be played into the GLOBAL variable E 
(GLOBAL so only one Vector.play() needs to be executed), 
which is multiplied by the RANGE variable d ("distance 
between the local node and the zero potential plane") 
and then by -1.

That product, called ex in this mechanism, is the extracellular 
potential at the local node, i.e. is used to drive local 
e_extracellular.

2. Reports the contribution of local i_membrane to the 
total signal that would be picked up by an extracellular 
recording electrode.  This is computed as the product of 
rx (transfer resistance between the local node and the 
monopolar recording electrode), i_membrane (called im in 
this mechanism), and the surface area of the local segment, 
and is reported as er.  The total extracellularly recorded 
potential is the sum of all er_xtrau over all segments in 
all sections, and is to be computed at the hoc level, 
e.g. with code like

func fieldrec() { local sum
  sum = 0
  forall {
    if (ismembrane("xtrau")) {
      for (x,0) sum += er_xtrau(x)
    }
  }
  return sum
}

Bipolar recording, i.e. recording the difference in potential 
between two extracellular electrodes, can be achieved with no 
change to either this NMODL code or fieldrec(); the values of 
rx will reflect the difference between the potentials at the 
recording electrodes caused by the local membrane current, so 
some rx will be negative and others positive.

3. Enables local storage of xyz coordinates interpolated from 
the pt3d data.  These coordinates are used by hoc code that 
computes the transfer resistance that couples the membrane 
to extracellular stimulating and recording electrodes.

4. Enables local storage of the distance d of internal 
nodes from the zero potential plane.

Prior to NEURON 5.5, the SOLVE statement in the BREAKPOINT block 
used METHOD cvode_t so that the adaptive integrators wouldn't miss 
the stimulus.  Otherwise, the BREAKPOINT block would have been called 
_after_ the integration step, rather than from within cvodes/ida, 
causing this mechanism to fail to deliver a stimulus current 
when the adaptive integrator is used.

With NEURON 5.5 and later, this mechanism abandons the BREAKPOINT 
block and uses the two new blocks BEFORE BREAKPOINT and  
AFTER BREAKPOINT, like this--

BEFORE BREAKPOINT { : before each cy' = f(y,t) setup
  ex = is*rx*(1e6)
}
AFTER SOLVE { : after each solution step
  er = (10)*rx*im*area
}

This ensures that the stimulus potential is computed prior to the 
solution step, and that the recorded potential is computed after.
ENDCOMMENT

NEURON {
	SUFFIX xtrau
	RANGE rx, er, d
	RANGE x, y, z
	GLOBAL E
	POINTER im, ex
}

PARAMETER {
	: default transfer resistance between recording electrode and node
	rx = 1 (megohm) : mV/nA
	x = 0 (1) : spatial coords
	y = 0 (1)
	z = 0 (1)
}

ASSIGNED {
	v (millivolts)
:	is (milliamp)
	E (volt/m) : field strength
	d (micron) : distance of node from zero potential plane
	ex (millivolts)
	im (milliamp/cm2)
	er (microvolts)
	area (micron2)
}

INITIAL {
:	ex = is*rx*(1e6)
	ex = -E*d*(1e-3)
	er = (10)*rx*im*area
: this demonstrates that area is known
: UNITSOFF
: printf("area = %f\n", area)
: UNITSON
}

: Use BREAKPOINT for NEURON 5.4 and earlier
: BREAKPOINT {
:	SOLVE f METHOD cvode_t
: }
:
: PROCEDURE f() {
:	: 1 mA * 1 megohm is 1000 volts
:	: but ex is in mV
:	ex = E*d*(1e-3)
:	er = (10)*rx*im*area
: }

: With NEURON 5.5 and later, abandon the BREAKPOINT block and PROCEDURE f(),
: and instead use BEFORE BREAKPOINT and AFTER BREAKPOINT

BEFORE BREAKPOINT { : before each cy' = f(y,t) setup
:  ex = is*rx*(1e6)
  ex = -E*d*(1e-3)
}
AFTER SOLVE { : after each solution step
  er = (10)*rx*im*area
}
