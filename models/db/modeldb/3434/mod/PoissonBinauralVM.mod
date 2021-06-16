TITLE PoissonBinauralVM.mod  Binaural Inhomogenous Poisson: 
 : binaural von Mises probability distributions

COMMENT
 Written by Jonathan Z Simon, jzsimon@isr.umd.edu
 A binaual inhomogenous poisson process
 Can be tailored to use many different distributions as stimuli.
 Stimulus Specific sections highlighted for makers of new stimuli 
 :------------------STIMULUS SPECIFIC Description --------------------------
  Uses 2 von Mises distributions (see e.g. Fisher, NI, Statistical Analysis
  of Circular Data), differing only in phase between ipsi and contra.
 :--------------------------------------------------------------------------
ENDCOMMENT

NEURON {

 :------------------STIMULUS SPECIFIC POINT_PROCESS name -------------------
   POINT_PROCESS PoissonBinauralVM
 :--------------------------------------------------------------------------

 RANGE 

 :------------------STIMULUS SPECIFIC RANGE variables ----------------------
  : PARAMETERs
   vs, freq, stimPhaseIpsi, stimPhaseContra, stimRate,
   genericParam1, genericParam2,
  : ASSIGNEDs
   angfreq, stimProb, kappa, norm, stim,
 :--------------------------------------------------------------------------
 
  poissDt, tickSize, tickVal, tickOnce : Stimulus independent
 POINTER probIpsi, probContra
}

UNITS {

 :------------------STIMULUS SPECIFIC UNITS --------------------------------
  (kHz) = (kiloHz)
 :--------------------------------------------------------------------------

 PI = (pi) (1)
 (tick) = (1)
 (prob) = (1)
 (angle) = (1) :should be radians but those are defined as .5 / pi
}

PARAMETER {

 :------------------STIMULUS SPECIFIC PARAMETERs ---------------------------
  vs = 0.43 (1)               <0,0.999> :otherwise too peaked even for exp()
  freq = 1000 (Hz)
  stimPhaseIpsi = 0 (angle)
  stimPhaseContra = 0 (angle)
  stimRate = 3.46 (/ms)       <1e-9,1e9>
  genericParam1 = 0 : unused
  genericParam2 = 0 : unused
 :--------------------------------------------------------------------------

 poissDt = 0.025 (ms)        <1e-9,1e9>
 tickSize = 1 (tick)         <1e-9,1e9>
}

ASSIGNED {
 :------------------STIMULUS SPECIFIC ASSIGNEDs ----------------------------
  angfreq (/ms)
  stimProb (prob)
  kappa (1)
  norm (1)
 :--------------------------------------------------------------------------

 tickVal (tick)
 tickOnce (tick)
 stim[2]  (prob) : [0] is Ipsi & [1] is contra
 probIpsi (prob)
 probContra (prob)

}

INITIAL {

 :------------------STIMULUS SPECIFIC INITIALizations ----------------------
  angfreq = 2*PI*freq*(0.001) : (0.001) from Hz & ms
  stimProb = stimRate*poissDt
  if (vs < 0) {vs = 0} 
  if (vs >= 0.999) {vs = 0.999} :otherwise too peaked even for exp()
  kappa =  A1Inv(vs)
  norm = 1/(2*PI*IO(kappa))
 :--------------------------------------------------------------------------

 tickVal = 0
 tickOnce = 0

 stim[0] = calcStimIpsi()
 stim[1] = calcStimContra()
 if (pointersValid()) {
  probIpsi = stimProb*stim[0]
  probContra = stimProb*stim[1]  
 }
}

:------------------STIMULUS SPECIFIC Stimulus Functions -------------------
 
 FUNCTION calcStimIpsi() {

  :von Mises probablility distribution (phase = 2 pi freq t + phiIpsi)

  calcStimIpsi = norm*exp(kappa*sin(angfreq*t+stimPhaseIpsi))

 }
 
 FUNCTION calcStimContra() {

  :von Mises probablility distribution (phase = 2 pi freq t + phiContra)

  calcStimContra = norm*exp(kappa*sin(angfreq*t+stimPhaseContra))

 }
 
:--------------------------------------------------------------------------

:------------------STIMULUS SPECIFIC Misc. Functions & Procedures----------
 FUNCTION IO(x) { 
 : Bessel Function I_O(x)
 : approximation from Fisher, NI, Statistical Analysis of Circular Data
  LOCAL w, w2, w4, w6, w8 : local vars for computational speed only
  w = x/3.75
  w2 = w*w
  w4 = w2*w2
  w6 = w4*w2
  w8 = w4*w4
  if (w < 1) {
   IO = (1 + 3.5156229*w2 + 3.0899424*w4 + 1.2067492*w6 + 0.2659732*w8 + 
             0.0360768*w8*w2 + 0.0045813*w8*w4)
  } else {
   IO = (0.39894228 + 0.01328592/w + 0.00225319/w2 - 0.00157565/(w2*w) + 
         0.00916281/w4 - 0.02057706/(w4*w) + 0.02635537/w6 - 
         0.01647633/(w6*w) + 0.00392377/w8)*exp(x)/sqrt(x)
  }
 }
 
 FUNCTION A1Inv(x) { 
 : Bessel Function Ratio: A_1(x) = I_1(x)/I_O(x), and use Inverse
 : approximation from Fisher, NI, Statistical Analysis of Circular Data
  if (x < 0.53) {
   A1Inv = (2*x + x*x*x + 5/6*x*x*x*x*x)
  } else { 
   if (x < 0.85) {
    A1Inv = (-0.4 + 1.39*x + 0.43/(1-x))
   } else {
    A1Inv = (1/(x*x*x - 4*x*x + 3*x))
   }
  }
 }

:--------------------------------------------------------------------------

: Inelegant but necessary to get poisson "clock" started.
: Must be called once after finitialize() but before run
PROCEDURE prerun() {
 tickVal  = tickSize
 tickOnce = tickSize
}

NET_RECEIVE (weight) {
 if(tickVal) {
  tickVal = 0
 } else {
  tickVal = tickSize
  stim[0] = calcStimIpsi()
  stim[1] = calcStimContra()
  probIpsi   = stimProb*stim[0]
  probContra = stimProb*stim[1]
 }
 net_send(poissDt,1)
}

FUNCTION pointersValid() {
 if ((!nrn_pointing(probIpsi))||(!nrn_pointing(probContra))) {
  printf(
   "Poisson Process: probIpsi and probContra must be set using setpointer!\n")
  pointersValid = 0
 } else {
  pointersValid = 1
 }
}

