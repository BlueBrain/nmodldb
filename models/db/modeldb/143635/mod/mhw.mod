COMMENT

Measures peak depol and calculates spike half width from the times at
which v crosses a (depolarized) threshold.  Threshold may be specified
by the user, or determined in the previous run.

Ted Carnevale

20110616 Modified by Tom Morse to follow the prescriptions of Amanda
Casale's matlab code threshdetect_neuronmod.m The notation was changed
to match the matlab variable names.  Also used method from Arnd Roth's
peak.mod from modeldb accession #135838 for dvdt.

20110522 Modified by Tom Morse to add measurements of 10-90% rise and
fall times as suggested by Amanda Casale.  Also the minimum voltage
used for calculations was changed from the initialization voltage to
an actual minimum (if run in mode 1,2)

USAGE EXAMPLES

//////////////////////////
// User-specified threshold
forall insert mhw
forall for (x,0) Halfheight_mhw(x) = THRESH // must assign value everywhere
mode_mhw = 0 // determine half width from fixed threshold
run()
printf(" base \t peak \t Halfheight \thalf width\n")
printf("%6.2f \t%6.2f \t%6.2f \t%6.2f\n", \
       Vthresh_mhw(0.5), Peak_mhw(0.5), Halfheight_mhw(0.5), Width_mhw(0.5))
//////////////////////////

//////////////////////////
// Dynamically-determined threshold
// run two simulations, first time with parameter mode_mhw = 1
//   and second time with mode_mhw = 2
// At end of first run, tmax and Peak will equal time and value of peak depol, and
// Vthresh will be set to the v where the membrane v's derivative dvdt exceeded
//    slope_thresh, a parameter set by the user prior to calling.
// At end of second run, Halfheight will be threshold for measurement of spike half width,
//    Risetime and Falltime will be the 90% rise time and fall time respectively,
//    t0 and t1 will be threshold crossing times, and Width will be spike half width
forall insert mhw
mode_mhw = 1
run() // find local Peak and associated tmax, thresh
mode_mhw = 2
run() // find spike Width
printf(" base \t peak \t Halfheight \thalf width\n")
printf("%6.2f \t%6.2f \t%6.2f \t%6.2f\n", \
       Vthresh_mhw(0.5), Peak_mhw(0.5), Halfheight_mhw(0.5), Width_mhw(0.5))
//////////////////////////

Width will be -1 if there is no max, or if simulation ends before v falls below Halfheight

Be cautious when using with adaptive integration--if the integrator uses long dt, 
t0 or t1 may be missed by a wide margin.
ENDCOMMENT

NEURON {
  SUFFIX mhw
  : mode values--
  : fixed threshold--0 use user-specified Halfheight
  : dynamic threshold--1 measure Peak, 2 calc Halfheight and measure halfwidth
  GLOBAL mode, slope_thresh
  RANGE Vthresh, Peak, tmax, Amp, dvdt : dvdt is actually "dV/dt"
  RANGE Halfheight, t0, t1, Width : Width is the half width
  RANGE Risetime, Falltime, t_rt_0, t_rt_1, t_ft_0, t_ft_1, thresh10, thresh90
}

UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
  (mM) = (milli/liter)
: note that ms is predefined in NEURON as 0.001 second
 }

PARAMETER {
  : mode values--
  : fixed threshold--0 use user-specified Halfheight
  : dynamic threshold--1 measure Peak, 2 calc Halfheight and measure halfwidth
  mode = 0 (1) : default, 0, is fixed (user-specified) threshold; (1,2) used in two runs
  slope_thresh = 20 (1) : 20 (mV/ms) arbitrary suggested value used to trigger the finding
    : of the approximate "threshold" of the neuron
}

ASSIGNED {
  v (mV)     : local v
  Vthresh (mV) : initial local v
  VthreshSet (1) : set to true once Vthresh is found in mode 1
  Peak (mV)  : max local v during previous run
  tmax (ms)  : time at which Peak occurred
  Halfheight (mV) : (Vthresh + Peak)/2
  t0 (ms)    : time in rising phase of spike when v first > Halfheight
  t1 (ms)    : time in falling phase of spike when v first < Halfheight
  Amp (mV)   : amplitude of AP from threshold to peak
  Width (ms)    : t1-t0
  findwhich (1) : assigned by program: 0 to find t0, 1 to find t1, 2 to find neither
  Risetime (ms)    : rise time: delta t between 10% and 90% of spike during rise to peak
  Falltime (ms)    : fall time: delta t between 90% and 10% of spike during fall from peak
  thresh10 (mV) : the 10% peak voltage threshold
  thresh90 (mV) : the 90% peak voltage threshold
  t_rt_0 (ms) : time of crossing thresh10 while rising
  t_rt_1 (ms) : time of crossing thresh90 while rising
  t_ft_0 (ms) : time of crossing thresh90 while falling
  t_ft_1 (ms) : time of crossing thresh10 while falling
                    : below help find t_rt_0, t_rt_1, t_ft_0, t_ft_1
  rise10set (1) : these all start at 0 (false) and then turn true (1) after set
  rise90set (1) : note the 1's here just mean these are unitless, in fact, boolean.
  fall90set (1)
  fall10set (1)
  
  v2 (mV)
  v3 (mV)
  t2 (mV)
  t3 (mV)
  dvdt (mV/ms)
}

INITIAL {
  if (mode==1) { : measure peak v then calc Halfheight
: printf("Finding Peak\n")
    Vthresh = -100 (mV) : hopefully will be overwritten with actual threshold
    Peak = v
    tmax = -1 (ms) : nonsense values for tmax, t0, t1, Width
    Halfheight = v
    t0 = -1 (ms)
    t1 = -1 (ms)
    Width = -1 (ms)
    Risetime = -1 (ms)
    Falltime = -1 (ms)
    thresh10 = v :these are reset to threshold at start of mode==2
    thresh90 = v
    t_rt_0 = -1 (ms)
    t_rt_1 = -1 (ms)
    t_ft_0 = -1 (ms)
    t_ft_1 = -1 (ms)
    : variables for calculating dvdt follow
	v2=0 (mV)
	v3=0 (mV)
	t2=0 (ms)
	t3=0 (ms)
	dvdt=0 (mV/ms)
    VthreshSet=0 (1) : gets set once Vthresh is found
	} else if (mode==2) { : calc Halfheight from Vthresh and Peak in order to determine halfwidth
: printf("Determining depolarization halfwidth\n")
    Amp = Peak - Vthresh
    Halfheight = Vthresh + 0.5 * Amp : there are a lot of ways to write this however this way is consistent with below
    findwhich = 0 : 0 to find t0, 1 to find t1
    thresh10 = Vthresh + 0.1 * Amp
    thresh90 = Vthresh + 0.9 * Amp
  } else if (mode==0) {
    Vthresh = v
    Peak = v
    tmax = -1 (ms) : nonsense values for tmax, t0, t1, Width, Risetime, Falltime, threshX, t_Xt_Y
    t0 = -1 (ms)
    t1 = -1 (ms)
    Width = -1 (ms)
    Risetime = -1 (ms)
    Falltime = -1 (ms)
:    thresh10 = hoc code user assigned
:    thresh90 = hoc code user assigned
    t_rt_0 = -1 (ms)
    t_rt_1 = -1 (ms)
    t_ft_0 = -1 (ms)
    t_ft_1 = -1 (ms)
    findwhich = 0 : 0 to find t0, 1 to find t1
  }
  rise10set = 0 (1) : in all modes can set false because rise fall time settings in
  rise90set = 0 (1) : findx only called in modes 0, 2 and not mode 1 where threshX unknown
  fall90set = 0 (1)
  fall10set = 0 (1)
 
  }

PROCEDURE findmax() {
  if (v>Peak) {
    Peak = v
    tmax = t
  }
  if (!VthreshSet && t > 2) { : t>2 avoids accidental current clamp signal
    if (dvdt >= slope_thresh) {
      Vthresh = v
      VthreshSet = 1
     }
  }
}

: find threshold crossings
PROCEDURE findx() {
: for half width
  if (findwhich==0) {
    if (v > Halfheight) {
      t0 = t
      findwhich = 1
    }
  } else if (findwhich==1) {
    if (v < Halfheight) { 
      t1 = t
      Width = t1-t0
      findwhich = 2 : stop looking already
    }
  }
: for fall times
  if (rise90set) { : only test after the rise times are set
    if (!fall90set) {
      if (v < thresh90) {
        t_ft_0 = t
        fall90set = 1
      }
    }
    if (!fall10set) {
      if (v < thresh10) {
        t_ft_1 = t
        fall10set = 1
        Falltime = t_ft_1 - t_ft_0
      }
    }
  }

: for rise times
  if (!rise10set) {
    if (v > thresh10) {
      t_rt_0 = t
      rise10set = 1
    }
  }
  if (!rise90set) {
    if (v > thresh90) {
      t_rt_1 = t
      rise90set = 1
      Risetime = t_rt_1 - t_rt_0
    }
  }
}

BREAKPOINT {
        SOLVE check METHOD after_cvode
}
PROCEDURE check() {
	v2 = v3
	v3 = v
	t2 = t3
	t3 = t
	if (t3 - t2 > 0) {
		dvdt = (v3 - v2)/(t3 - t2)
	}
}
AFTER SOLVE { : works as well, executed half as many times
  if (mode==1) { : measure peak v (Peak) and Vthresh
    findmax()
  } else if (mode==2) {
    findx()
  } else if (mode==0) {
    findmax() : might as well, even if we don't use it to find threshold
    findx()
  }
}
