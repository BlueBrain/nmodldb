: $Id: presyn.inc,v 1.1 2005/04/04 23:35:42 hines Exp hines $

NEURON {
  RANGE maxsyn, nsyn                         : scalars
  RANGE link                              : PreL structure
}

PARAMETER {
  nsyn = 0
  link = 0       : presyn struct points to list of pointers to postsyn cells
}

ASSIGNED {
  maxsyn        : max and counter for the list of pointers
}

INCLUDE "snshead.inc"

CONSTRUCTOR {
  VERBATIM {
    /* link allows postsyn cell to manipulate presyn list */
    link = (double)((unsigned long)hoc_Ecalloc(1, sizeof(PreL))); hoc_malchk();
    PRECAST->link2  = link;  /* this can be queried to protect against error */
    PRECAST->nnpre = &nsyn;
    PRECAST->mxpre = &maxsyn;
    nsyn = 0.;
    maxsyn = 0.;
    if (ifarg(2)) { 
      PRECAST->cpre = (int)*getarg(2);
    } else { PRECAST->cpre = -1; }
 } 
 ENDVERBATIM
}

DESTRUCTOR {
  VERBATIM { 
      nsyn = maxsyn = 0;
      free(PRECAST);
      link = 0;
  }
  ENDVERBATIM
}

PROCEDURE preinit () {}

: globally used presynaptically to add to queue (was called pushqueu())
PROCEDURE newspike() {
VERBATIM { 
  int ii, jj, nspk, qend, head, ind, swapi;
  double time, swapt;
  SynS *psyn;		/* pointer to the presynaptic structure */
  QptR *qpr;
  QueU *qp;		/* pointer to the head position in postsyn queu */

  for (ii=0; ii<nsyn;ii++) { 

    /* psyn is the pointer back to the postsyn info array (struct SynS) */
    if ((psyn = PRECAST->plst[ii])==(SynS *)NULL ||     /* skip empties */
	(t < psyn->spkt && psyn->chainlen == 1)) {   /* skip if still on last one */
      continue; }
    /* are we chaining synapses? : this means that we will treat additional spks as */
    /* if they were coming from another synapse */
    if (psyn->chainlen > 1) {
      psyn->chainptr++;  /* augment the pointer */
      if (psyn->chainptr >= psyn->chainlen) { /* have run out of room */
	if (t < psyn->spkt) {  /* can't reset yet so just drop this spk */
	  continue;
	} else {               /* the whole chain is completely done so can reset */
	  psyn->chainptr = 0; }
      }
      ind = psyn->index + psyn->chainptr; /* add index of next one on chain */
    } else {
      ind = psyn->index;
    }

    /* qpr = the post-syn queu pointer */
    qpr = psyn->qpt;  
    /* qp = pointer to postsyn queu */
    qp = qpr->qq;
    nspk = (int)(*(qpr->nspk));
    time = t + psyn->del;
    psyn->spkt = t + qpr->dead; /* when it will all end */

    /* basic queu code taken from queu.inc (viz.) but now must make sure */
    /* that we don't go backwards in time (different delays) */
    /* should have been cleared if we are to write here */
    if (qp[nspk].time == 1e20) {
      if (qpr->qterm < time) {	/* queu marching forward in time: AOK */
	/* push the two queu items: a time and an index */
	(qp[nspk]).time = time; /* delay */
	(qp[nspk]).index = ind; /* tell postsyn which presyn spiked */
	qpr->qterm = time;	/* keep track of largest time put on */
      } else {			/* move backward to insert in right place */
	/* actually a heap rather than a queue since maintains order */
	qend = (int)(*(qpr->qln)) - 1; /* max ind for queu */
	head = (int)(*(qpr->head)) - 1; /* one in front of the head */
	if (head == -1) { head = qend; } /* wrap around */
	/* search backward */
	for (jj=((nspk==0)?qend:nspk-1);
	     time < qp[jj].time && jj != head;
	     jj=((jj==0)?qend:jj-1)); /* skip over from end */
	if (jj==nspk) { 
	  printf("%g %d %d ",qpr->qterm,ii,jj);
	  hoc_execerror("ERROR: newspike search failed.\n",0); }
	/* move things forward */
	for (jj=((jj==qend)?0:jj+1); /* start one up */
	     jj != nspk;	/* till back where we were */
	     jj=((jj==qend)?0:jj+1)) {
	  swapt = qp[jj].time;	/* store temporarily */
	  swapi = qp[jj].index;
	  qp[jj].time = time;	/* put in the new one */
	  qp[jj].index = ind;
	  time = swapt;		/* swap them */
	  ind = swapi;
	}
	qp[jj].time = time;
	qp[jj].index = ind;
      }

      /* move the tail forward and reset if reached the end */
      (*(qpr->nspk))++;
      if ((*(qpr->nspk)) == *(qpr->qln)) { (*(qpr->nspk)) = 0.; }

    } else {			/* (qp->time != 1e20) */
      printf("%x %g %g ",qpr->qq,*(qpr->nspk),*(qpr->qln));
      hoc_execerror("Error: queue full.\n",0);
    }
  }
  }
ENDVERBATIM
}

: eliminates empty pointers in presyn array, moving everything else down
PROCEDURE clean() {
  VERBATIM { 
    int ii,jj,newn;
    SynS** pl;
    pl = PRECAST->plst;

    for(ii=0,jj=0,newn=nsyn; ii<nsyn; ii++,jj++) {
      while (pl[ii] == (SynS *)NULL && ii<nsyn) {
	ii++; newn--; 
      }
      if (ii != jj && ii<nsyn) {
	pl[jj] = pl[ii];
	pl[jj]->back = &(pl[jj]);
      }
    }
    nsyn = newn;
  }
  ENDVERBATIM
}

: prints out some selected information from the presyn structure, synlist
FUNCTION check() {
  VERBATIM { 
    int ii;
    for (ii=0,_lcheck=0; ii<nsyn && _lcheck!=-1.; ii++) {
      if (PRECAST->plst[ii] != (SynS *)NULL) {
	_lcheck++;  /* count the active pointers */
	if ((&(PRECAST->plst[ii])) != (PRECAST->plst[ii]->back)) {
	  printf("****************************************************************\n");
	  printf("ERROR:: %2d: Index:%3d,C%09d,Pr%09d,Po%09d, Delay:%6g (%x->%x->%x)\n",
		 ii,
		 PRECAST->plst[(int)ii]->index,
		 PRECAST->plst[(int)ii]->ucode,
		 PRECAST->cpre,
		 PRECAST->plst[(int)ii]->qpt->cpost,
		 PRECAST->plst[(int)ii]->del,
		 &(PRECAST->plst[(int)ii]),
		 PRECAST->plst[(int)ii],
		 (PRECAST->plst[(int)ii]->back));
	  printf("****************************************************************\n");
	  _lcheck = -1.0;
	}
      }
    } 
  }
  ENDVERBATIM
}

: return ucode
FUNCTION code() {
  VERBATIM { int ii;
    if (ifarg(1)) {
      ii = (int)*getarg(1);
      if (ii < 0) { ii = nsyn+ii; }
      if (ii >= nsyn || ii < 0) {hoc_execerror("array index out of bounds", 0);}
      _lcode = ((PRECAST->plst[ii] == (SynS *)NULL)?-2:PRECAST->plst[ii]->ucode);
    } else {
      for (ii=0;ii<nsyn;ii++) {
	if (PRECAST->plst[ii] == (SynS *)NULL) {
	  printf("%d EMPTY.\n",ii);
	} else { 
	  printf("%d,%09d  ",ii,PRECAST->plst[ii]->ucode);
        }
      }
      _lcode = 1.0;
    }
  } 
  ENDVERBATIM
}

: return delay
FUNCTION delay() {
  VERBATIM { int ii;
    if (ifarg(1)) {
      ii = (int)*getarg(1);
      if (ii < 0) { ii = nsyn+ii; }
      if (ii >= nsyn || ii < 0) {hoc_execerror("array index out of bounds", 0);}
      _ldelay = ((PRECAST->plst[ii] == (SynS *)NULL)?-2:PRECAST->plst[ii]->del);
    } else {
      for (ii=0;ii<nsyn;ii++) {
	if (PRECAST->plst[ii] == (SynS *)NULL) {
	  printf("%d EMPTY.\n",ii);
	} else { 
	  printf("%d,%g  ",ii,PRECAST->plst[ii]->del);
        }
      }
      _ldelay = 1.0;
    }
  } 
  ENDVERBATIM
}

: return or set cpre code
FUNCTION pre() {
  VERBATIM {
    if (ifarg(2)) { PRECAST->cpre = (int)*getarg(2); }
    _lpre = PRECAST->cpre;
  } 
  ENDVERBATIM
}

: return cpost code
: 0 arg: prints out information from the postsyn structure
: 1 arg = -1 print out the pointer structure
: plst[i] remains only indication of eliminated pointer
FUNCTION post() {
  VERBATIM { int xx,ii;
    xx = -2;
    if (ifarg(1)) { xx = (int)*getarg(1); }
    if (xx >= 0) {
      if (xx >= nsyn) {hoc_execerror("array index out of bounds", 0);}
      _lpost = ((PRECAST->plst[xx] == (SynS *)NULL)?-2:PRECAST->plst[xx]->qpt->cpost);
    } else {
      for (ii=0;ii<nsyn;ii++) {
	if (PRECAST->plst[ii] == (SynS *)NULL) {
	  printf("%2d EMPTY.\n",ii);
	} else { 
	  printf("%1d (%1d) Del:%-4g Cd:%09d Pre:%09d Post:%09d\n",
		 ii,
		 PRECAST->plst[ii]->index,
		 PRECAST->plst[ii]->del,
		 PRECAST->plst[ii]->ucode,
		 PRECAST->cpre,
		 PRECAST->plst[ii]->qpt->cpost);
	  if (xx == -1 || (&(PRECAST->plst[ii])) != (PRECAST->plst[ii]->back)) {
	    printf("\t\t%s: %x -> %x\n",
		   (((&(PRECAST->plst[ii]))==(PRECAST->plst[ii]->back))?
		    "OK":"POINTER MISMATCH"),
		   &(PRECAST->plst[ii]),
		   PRECAST->plst[ii]);
	  }
	}
      }
      _lpost = 1.0;
    }
  } 
  ENDVERBATIM
}

: for debugging - prints out queu info from the pre side
PROCEDURE prq(x) {
  VERBATIM { int ii;
    SynS *psyn;		
    QptR *qpr;
    psyn = PRECAST->plst[(int)_lx];
    qpr = psyn->qpt;  
    printf("%x %g %g\n",qpr->qq,*(qpr->nspk),*(qpr->qln));
    for (ii=0;ii<*(qpr->qln);ii++) {
      printf("%8g%8d\n",
	     ((qpr->qq)[ii]).time,
	     ((qpr->qq)[ii]).index);
   }
  } 
  ENDVERBATIM
}

