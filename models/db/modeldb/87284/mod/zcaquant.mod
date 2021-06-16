COMMENT
caquant.mod
based on mhw.mod, which measures peak depol and calculates spike half width
  from the times at which v crosses a (depolarized) threshold.

caquant measures these:
iinit   ica at t==0, in mA/cm2
imax    most negative (i.e. largest inward) ica, in mA/cm2
timax   time at which the most negative ica occurs
hwi     half width of the first ica transient that occurs in a simulation
          and has most negative value (i.e. dips below) < imax + irest
cinit   cai at t==0
cmax    largest cai
tcmax   time at which the largest cai occurs

caquant also calculates these values:
qapprox = -area*(imax-iinit)*hwi
        approximate "extra" charge that enters during the ica transient
        (assumes waveform is more or less triangular)
        ("extra" charge means the difference between the amount that enters 
        during the ica transient and the amount that would have entered 
        if the transient hadn't happened)
cmaxp   predicted cmax calculated as cinit - svr*(imax-iinit)*hwi/2/FARADAY
        where svr is a RANGE variable (a PARAMETER declared to be RANGE) 
        equal to the surface/volume ratio.  svr must be specified by the 
        user during setup, because it depends on the geometry assumed by 
        the accumulation mechanism that READs ica and WRITEs cai.

caquant operates in two modes, which are controlled by the PARAMETER (a global)
called mode.
When mode is 1, caquant determines the following in the course of a simulation:
iinit, imax, timax, cinit, cmax, and tcmax
When mode is 2, it does the following in the course of a simulation:
1.  uses the values of iinit and imax to set the threshold that it monitors 
    in order to determine hwi  
2.  uses the values of hwi, iinit, and imax to calculate qapprox,
and
3.  uses the values of cinit, svr, imax, iinit, and hwi to calculate cmaxp.

USAGE EXAMPLE

/*
Run two simulations, first time with parameter mode_caquant = 1
and second time with mode_caquant = 2

At end of first run,
this  will equal  this
iinit             ica at t==0, in mA/cm2
imax              most negative (i.e. largest inward) ica, in mA/cm2
timax             time at which the most negative ica occurs
cinit             cai at t==0
cmax              largest cai
tcmax             time at which the largest cai occurs

At end of second run,
this  will equal  this
t0i and t1i        the times at which ica crosses the "threshold" iinit+imax/2
hwi              half width of the first ica transient
qapprox          approximate extra charge that enters during the ica transient
cmaxp            predicted cmax
*/

forall {
  insert caquant
  for (x, 0) svr_caquant(x) = some function that specifies the surface/volume ratio
    of the current segment
}
mode_caquant = 1
run()
mode_mhw = 2
run()

hwi will be -1 if there is no max, or if simulation ends before ica crosses
  the automatically-set threshold twice

Be cautious when using with adaptive integration--if the integrator uses long dt, 
t0i or t1i may be missed by a wide margin.
ENDCOMMENT

NEURON {
  SUFFIX caquant
  USEION ca READ ica, cai
  : mode values--
  : 1 measure imax, 2 calc ihalf and measure halfwidth
  GLOBAL mode
  RANGE iinit, imax, timax
  RANGE ihalf, t0i, t1i, hwi
  RANGE cinit, cmax, tcmax

  RANGE svr : local surface/volume ratio
  RANGE qapprox, cmaxp

  RANGE vinit, vmax, tvmax
  RANGE vhalf, t0v, t1v, hwv
  RANGE vxt
}

UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
  (mM) = (milli/liter)
  F = (faraday) (coulombs)
}

PARAMETER {
  : mode values--
  : 1 measure imax, 2 calc ihalf and measure halfwidth
  mode = 1 (1) : default is measure imax
}

ASSIGNED {
  ica (mA/cm2)
  iinit (mA/cm2) : initial ica
  imax (mA/cm2)  : max ica during previous run
  timax (ms)  : time at which imax occurred
  ihalf (mA/cm2) : (iinit + imax)/2
  t0i (ms)    : time in rising phase of cai transient when cai first < ihalf
  t1i (ms)    : time in falling phase of cai when cai first > ihalf
  hwi (ms)    : t1i-t0i
  findwhichi (1) : 0 to find t0i, 1 to find t1i, 2 to find neither

  cai (mM)
  cinit (mM) : initial cai
  cmax (mM)  : max cai during previous run
  tcmax (ms)  : time at which imax occurred

  svr (/micron)  : surface/volume ratio, which depends on the ca accumulation mechanism
            : for cacum.mod it is depth_cacum

  area (micron2) : segment surface area
  qapprox (picocoulomb)
  cmaxp (mM)

  v (mV)
  vinit (mV) : initial v
  vmax (mV)  : max v during previous run
  tvmax (ms)  : time at which vmax occurred
  vhalf (mV) : (vinit + vmax)/2
  t0v (ms)    : time in rising phase of v transient when v first > vhalf
  t1v (ms)    : time in falling phase of v when v first < vhalf
  hwv (ms)    : t1v-t0v
  findwhichv (1) : 0 to find t0v, 1 to find t1v, 2 to find neither
  vxt (ms mV) : (vmax - vinit)*hwv
}

INITIAL {
  if (mode==1) { : measure peak ica
: printf("Finding imax\n")
    iinit = ica
: printf("iinit is %g\n", iinit)
    imax = ica
    timax = t
    ihalf = ica
    t0i = -1 (ms) : nonsense values for t0i, t1i, hwi
    t1i = -1 (ms)
    hwi = -1 (ms)

    cinit = cai
    cmax = cai
    tcmax = t

    vinit = v
    vmax = v
    tvmax = t
    vhalf = v
    t0v = -1 (ms) : nonsense values for t0v, t1v, hwv
    t1v = -1 (ms)
    hwv = -1 (ms)

  } else if (mode==2) { : calc ihalf from iinit and imax in order to determine halfwidth
: printf("Determining ica transient's halfwidth\n")
    ihalf = (iinit + imax)/2
: printf("iinit is %g and ihalf is %g\n", iinit, ihalf)
    findwhichi = 0 : 0 to find t0i, 1 to find t1i
    vhalf = (vinit + vmax)/2
    findwhichv = 0 : 0 to find t0v, 1 to find t1v
  }
}

: find ica threshold crossings
PROCEDURE findix() {
  if (findwhichi==0) {
    if (ica < ihalf) {
      t0i = t
      findwhichi = 1
    }
  } else if (findwhichi==1) {
    if (ica > ihalf) { 
      t1i = t
      hwi = t1i-t0i
      findwhichi = 2 : stop looking already
    }
  }
}

: find v threshold crossings
PROCEDURE findvx() {
  if (findwhichv == 0) {
    if (v > vhalf) {
      t0v = t
      findwhichv = 1
    }
  } else if (findwhichv == 1) {
    if (v < vhalf) { 
      t1v = t
      hwv = t1v-t0v
      findwhichv = 2 : stop looking already
    }
  }
}

COMMENT
: steady state m from car mechanism
FUNCTION minf(V (mV))(1) {
  minf = 1 / (1 + exp((v+14(mV))/(-6.7 (mV))))
}
ENDCOMMENT

: BREAKPOINT {
: a mechanism that calculate something on every time step
: needs a BREAKPOINT block--even if only an empty one--
: or else it will be treated as an ARTIFICIAL_CELL
BREAKPOINT { }

: AFTER SOLVE { : should work as well, executed half as many times
BEFORE STEP { : should work even with cvode
  if (mode==1) { : measure peak ica, cai, v, and pcar
    if (ica<imax) {
      imax = ica
      timax = t
    }
    if (cai>cmax) {
      cmax = cai
      tcmax = t
    }
    if (v>vmax) {
      vmax = v
      tvmax = t
    }
  } else if (mode==2) {
    findix()
    qapprox = -(0.01)*area*(imax-iinit)*hwi
    cmaxp = cinit - (10000)*svr*(imax-iinit)*hwi/2/F
    findvx()
    vxt = (vmax - vinit)*hwv
  }
}
