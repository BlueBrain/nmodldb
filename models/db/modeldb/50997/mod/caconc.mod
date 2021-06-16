COMMENT
Mechanism for varying internal [Ca++]
ENDCOMMENT
					       
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS caconc
	RANGE  Alpha_Max, Alpha_Delay, Alpha_tau, DC_Level, DC_Delay, DC_Off, Ramp_Max, Ramp_Delay, Ramp_Off, Slope_UP, Slope_DOWN, alpend
	GLOBAL caconc
}

UNITS {
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {
        Alpha_Max=0     (umho)
	Alpha_Delay=0 (ms)
	Alpha_tau=.1 (ms)
	e=0	(mV)
	v	(mV)
	caconc=0	(mM)
        DC_Level (mM)
        DC_Delay (ms)
        DC_Off (ms)
        Ramp_Max (mM)
        Ramp_Delay (ms)
        Ramp_Off (ms)
        Slope_UP (ms)
        Slope_DOWN (ms)
        alpend=0
}


BREAKPOINT {
       if (Ramp_Max > 0)
          { if (Ramp_Off > t) 
              {caconc = gramp(t) + (Alpha_Max * alpha( (t - Alpha_Delay)/Alpha_tau ))
               alpend = caconc }
            else
              {caconc = alpend + (Slope_DOWN * (t - Ramp_Off))}
          }
       else
          {caconc = gramp(t) + (Alpha_Max * alpha( (t - Alpha_Delay)/Alpha_tau ))}

       if (caconc < 0) { caconc = 0 }


}

FUNCTION alpha(x) {
	if (x < 0 || x > 10) {
		alpha = 0
	}else{
		alpha = x * exp(1 - x)
	}
}



FUNCTION gramp(x)
{
VERBATIM
double tramp, x, Dc;

x = _lx;

if (x >= DC_Delay & x < DC_Off)
   Dc = DC_Level;
else
   Dc = 0;

if (x < Ramp_Delay)
   tramp = Dc;
else
  {
  if (x < Ramp_Off)
   {
    tramp = Dc + (Slope_UP * (x - Ramp_Delay));
    if (tramp >= Ramp_Max)
       tramp = Ramp_Max;
   }
  else
   {
    tramp = Ramp_Max + (Slope_DOWN * (x - Ramp_Off));
    if (tramp <= Dc) 
       tramp = Dc;
   }
  }
return (tramp);
ENDVERBATIM
}
