TITLE PoissonTrigSyn.mod   Poisson triggerable alpha synapse

COMMENT
 Poisson triggerable alpha function synapse
 written by Jonathan Z Simon, jzsimon@isr.umd.edu
 with core code taken from Hines & Carnivale, "Expanding NEURON with NMODL"
ENDCOMMENT

COMMENT
 modified by Jonathan Simon, jzsimon@isr.umd.edu
 a synaptic current with alpha function conductance that works out to be
  i = g * (v - e)
 where
  g = 0 before the synapse is triggered at onset
  g = stim * (t - onset)/tau * exp(-(t - onset - tau)/tau) for t > onset
 this has the property that the maximum value is gmax and occurs at
  t =  onset + tau.
 The onset time is triggered by a pseudo poisson process.
ENDCOMMENT

NEURON {
 POINT_PROCESS PoissonTrigSyn
 RANGE tau, refraction, gmax, e, i, g, k, firing
 POINTER prob, randseed
 NONSPECIFIC_CURRENT i
}

UNITS {
 (nA) = (nanoamp)
 (mV) = (millivolt)
 (uS) = (microsiemens)
}

PARAMETER {
 tau=.1 (ms)
 refraction = 1 (ms)
 gmax=0  (uS)
 e=0 (mV)
 randseed
}

ASSIGNED {
 i (nA)
 g (uS)
 k (/ms)
 firing
 prob
 v (mV)
}

CONSTANT { exp1 = 2.7182818284590452354} :Euler constant for correct max peak

STATE { A (uS) G (uS) }

INITIAL {
 k = 1/tau
 A = 0
 G = 0
 firing = 0
 if (!nrn_pointing(prob)) {
  printf("PoissonTrigSyn: prob must be set using setpointer!\n")
 }
 if (!nrn_pointing(randseed)) {
  printf("PoissonTrigSyn: randseed must be set using setpointer!\n")
 }
}


BREAKPOINT {
 SOLVE state METHOD sparse
 i = G*(v - e)
 g = G
}

KINETIC state {
 ~ A <-> G (k, 0)
 ~ G ->    (k)
}

NET_RECEIVE (weight) {
 if (flag == 0) {  : event generated externally--fire if not firing now
  if ((prob > unitRand()) && (!firing)){
   state_discontinuity(A, A + gmax*exp1)
   firing = 1
   if (refraction >= 0) {
    net_send(refraction, firing) : send event to clear firing flag
   } else {
    firing = 0                   : or clear firing flag immediately
   }
  }
 } else {
  firing = 0  : event generated internally to clear firing flag
 }
}

FUNCTION getseedCD() {
 VERBATIM
#define RANDSEED_UINT (*((unsigned int*)(&randseed)))
#define GETSEEDCD (_lgetseedCD)
  GETSEEDCD = (double)RANDSEED_UINT; /* emulate getseedCD = randseed */
 ENDVERBATIM
}

PROCEDURE initseedCD(seed) {
 VERBATIM
#define RANDSEED_UINT (*((unsigned int*)(&randseed)))
#define SEED (_lseed)
  RANDSEED_UINT = (int)SEED;   /* emulate randseed = seed */
  /* should be (unsigned int) instead of (int) but MPW lacks the necessary
      library function to do the conversion. Not a big deal */ 
 ENDVERBATIM
}

FUNCTION unitRand() { : taken from by c's stdlib rand() function & man
 VERBATIM
#define RANDSEED_UINT (*((unsigned int*)(&randseed)))
#define UNITRAND (_lunitRand)
   /* emulate randseed = randseed * 1103515245 + 12345 */
   RANDSEED_UINT = RANDSEED_UINT * 1103515245 + 12345;
   /* emulate unitRand = ((randseed >> 16) & 0x7FFF)/(RAND_MAX+1.0); */
   UNITRAND = ((RANDSEED_UINT >> 16) & 0x7FFF)/(((unsigned int)(0x7FFF))+1.0);
 ENDVERBATIM
}
