COMMENT
Ribbon synapse version 1

activation follows Poisson rate (atau) time constant
alpha(t) = t/tau * exp( -t/tau )

ENDCOMMENT

VERBATIM
extern double du_dev0( );
extern double dexp_dev( );
ENDVERBATIM

NEURON {
  POINT_PROCESS 	ribbon1
  RANGE 		i, del, xp, yp, zp, dist, tau, sw
  NONSPECIFIC_CURRENT 	i
}

PARAMETER {
  del	= -1	(ms)	: if >= 0 will trigger synaptic event at t = del
  xp	= 0		: xyz location of synapse
  yp 	= 0
  zp 	= 0
  dist	= 0		: distance to recording section
  e 	= 0   	(mV)	: equib potential
  tau	= 1.3	(ms)	: time constant of alpha current
  atau	= 0	(ms)	: time constant of activation (default=0 => inactive )
  sw	= 1000	(pS)	: synaptic weight for self-event
  mp	= 1.0		: mp>=1 mean of geometric multiple release probability (geometric with prob 1/mp)
  idebug = 0
}

ASSIGNED {
  v (mV)
  i (nanoamp)
}

STATE {
  a (pS)
  g (pS)
}

UNITS {
 (mV) 	= (millivolt)
 (pS) 	= (picosiemens)
 PI	= (pi) (1)
}

INITIAL {
  if( del >= 0 ){
    net_send( del, 3 )
  }
  net_send( 1.0e-10, 1 )	: seems sometimes misses if I use zero
  g = 0
}

BREAKPOINT {
  SOLVE state METHOD sparse
  i = g*( v - e )*(1e-06)	: Convert to pS
}

KINETIC state {
  ~ a <-> g ( 1/tau, 0 )
  ~ g  ->   ( 1/tau )
}

NET_RECEIVE( weight (pS)) {
  LOCAL wait
  UNITSOFF
  if( idebug ) {
    printf( "Net_receive afhcsyn: t %g flag %g", t, flag )
    if( flag == 0 ){ printf( " weight %g", weight ) }
    printf( "\n" )
  }
  if( flag == 0 ){	: NetCon event (use weight from event)
    a = a + weight * exp(1) * igeom_dev(1/mp)
  }
  if( flag == 1 && atau > 0 ){	: self-event
    wait = dexp_dev( 1(ms)/atau )
    net_send( wait, 2 )
  }
  if( flag == 2 && atau > 0 ){	
    a = a + sw *exp(1)* igeom_dev(1/mp)
    wait = dexp_dev( 1/atau )
    net_send( wait, 2 )
  }
  if( flag == 3 ){	: self-event from del > 0
    a = a + sw *exp(1) * igeom_dev(1/mp)
  }
  UNITSON
}
