TITLE SpikeStats.mod  Spike Statistics

COMMENT
 Spike Statistics calculator which can be run by a NetCon
 written by Jonathan Z Simon, jzsimon@isr.umd.edu
ENDCOMMENT

NEURON {
 POINT_PROCESS SpikeStats
 RANGE binsnum, ignoreBefore, stimPeriod, poissDt, sinsum, cossum
 POINTER n, vs, phase, pst, sinphs, cosphs
}

UNITS {
 PI = (pi) (1)
}

CONSTANT {radToDeg = 57.295779513082322865} : = 180/pi

PARAMETER {
 binsnum = 80                 <0,1e9>
 ignoreBefore = 0 (ms)        <0,1e9>
 stimPeriod  = 1 (ms)         <0,1e9>
 poissDt = 0.0125 (ms)        <1e-9,1e9>
 : Should have pst, sinphs, cosphs here too, but pointers 
 : (and arrays) work inelegantly in PARAMETER block. 
 : For hoc usage consistency, n, vs, phase, are also made pointers.
 : All these must be assigned with setpointer.
}

ASSIGNED {
 sinsum
 cossum
 n      : pointer
 vs     : pointer
 phase  : pointer
 pst    : pointer/array
 sinphs : pointer/array
 cosphs : pointer/array
}

INITIAL {
 LOCAL phsi
 VERBATIM
  /* define C versions of variables to be used below */
#define PST    (&pst)
#define SINPHS (&sinphs)
#define COSPHS (&cosphs)
#define I_INT  ((int)_li)
#define PHSI   (_lphsi)
 ENDVERBATIM
 sinsum = 0
 cossum = 0
 if (pointersValid()) {
  n = 0
  vs = 0
  phase = -1e-9
  FROM i=0 TO binsnum-1 {
   phsi = 2*PI/binsnum*i
   VERBATIM
    PST[I_INT]    = 0.0;       /* emulate pst[i]    = 0         */
    SINPHS[I_INT] = sin(PHSI); /* emulate sinphs[i] = sin(phsi) */
    COSPHS[I_INT] = cos(PHSI); /* emulate cosphs[i] = cos(phsi) */
   ENDVERBATIM
  }
 }
}

NET_RECEIVE(change) {
 LOCAL i
 VERBATIM
  /* define C versions of variables to be used below */
#define PST    (&pst)
#define SINPHS (&sinphs)
#define COSPHS (&cosphs)
#define I_INT  ((int)_li)
 ENDVERBATIM
 if (t > ignoreBefore) {
  n = n + 1
  UNITSOFF : shouldn't need to turn unit checking off here, but do
   i = fmod(floor(fmod(t,stimPeriod)/poissDt),binsnum)
  UNITSON
  VERBATIM
   (PST[I_INT])++ ;         /* emulate pst[i] = pst[i] + 1         */
   sinsum += SINPHS[I_INT]; /* emulate sinsum = sinsum + sinphs[i] */
   cossum += COSPHS[I_INT]; /* emulate cossum = cossum + cosphs[i] */
  ENDVERBATIM
  vs = sqrt(sinsum*sinsum + cossum*cossum)/n
  phase = atan2(sinsum,cossum)*radToDeg
  if (phase < 0) { phase = phase + 360 }
 }
}

FUNCTION pointersValid() {
 if ((!nrn_pointing(n))||(!nrn_pointing(vs))||(!nrn_pointing(phase))||
     (!nrn_pointing(pst))||(!nrn_pointing(sinphs))||
     (!nrn_pointing(cosphs))) {
  printf("SpikeStats: n, vs, phase, pst, sinphs, and cosphs ")
  printf("must all be set using setpointer!\n")
  pointersValid = 0
 } else {
  pointersValid = 1
 }
}

