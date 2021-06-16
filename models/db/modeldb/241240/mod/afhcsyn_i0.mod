COMMENT
Alpha function Hair Cell synapse 

activation follows Poisson rate (atau) time constant
alpha(t) = t/tau * exp( -t/tau )

ENDCOMMENT

VERBATIM
extern double du_dev0( );
extern double dexp_dev( );
ENDVERBATIM

NEURON {
  POINT_PROCESS 	afhcsyn_i0
  RANGE 		i, tau, atau, sw, del
  NONSPECIFIC_CURRENT 	i
}

PARAMETER {
  e 	= 0   	(mV)	: equib potential
  tau	= 1.3	(ms)	: time constant of alpha current
  atau	= 0	(ms)	: time constant of activation (default=0 => inactive )
  sw	= 0.1	(pS)	: default synaptic weight
  del	= 0.1	(ms)	: default delay

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
  net_send( 0, 1 )
  g = 0
}

BREAKPOINT {
  SOLVE state METHOD sparse
  i = g*( v - e )*(1e-06)
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
  if( flag == 0 ){
    a = a + weight * exp(1)
  }
  if( flag == 1 && atau > 0 ){	
    wait = del                   
    net_send( wait, 2 )
  }
  if( flag == 2 && atau > 0 ){	
    wait = dexp_dev( 1(ms)/atau )
    net_send( wait, 3 )
  }
  if( flag == 3 && atau > 0 ){	
    a = a + sw *exp(1)
    wait = dexp_dev( 1/atau )
    net_send( wait, 1 )
  }
  UNITSON
}
