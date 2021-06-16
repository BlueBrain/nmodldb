: $Id: snsarr.inc,v 1.63 1998/10/07 15:06:12 billl Exp $

COMMENT
manage arrays used in sns.inc
ENDCOMMENT

NEURON {
  RANGE maxsyn, nsyn, qlen                   : scalars
  RANGE lpst, queu                            : SynS bzw. QueU structures
  RANGE begsyn, endsyn, newspk
  GLOBAL SPKSTORE, CHAINLEN
  RANGE qptr                                 : QptR structure
}

PARAMETER {
  SPKSTORE = 50                  : max number of spikes expected during max delay
  CHAINLEN = 1                  : how many extra entries in a chain
  lpst = 0                      : pointer to postsynaptic array
}

ASSIGNED {
  maxsyn        : max and counter for this array
  nsyn		
		
  queu          : the postsyn queu
  qlen          : its length (maxsyn*SPKSTORE)
  begsyn        : index into queu that tells time to start
  endsyn        : index into queu that tells time to end
  newspk        : index into queu to place spk time + delay (tail)

  qptr          : will point to the queu, its tail and its length
}

INCLUDE "snshead.inc"
:  QptR qptr;  /* governs access to the postsyn queue */
:  PreL prel;  /* governs access to the presyn list */

: takes 1-3 arguments: location, post code, maxsyn
CONSTRUCTOR {
  VERBATIM {
    qptr = (double)((unsigned long)ecalloc(1, (sizeof(QptR)))); 
    /* qptr allows presynaptic cell to easily manipulate the post queu */
    QPRCAST->qln = &qlen;
    QPRCAST->nspk = &newspk;
    QPRCAST->head = &begsyn;
    nsyn = maxsyn = qlen = lpst = queu = 0.;
    if (ifarg(2)) { 
      QPRCAST->cpost = (int)*getarg(2);
    } else { QPRCAST->cpost = -1; }
    if (ifarg(3)) { init_arrays(*getarg(3)); }
  }
  ENDVERBATIM
}

DESTRUCTOR {
  VERBATIM { int ii;
      free(QUECAST);  
      if (lpst != 0) {
        for (ii=0;ii<nsyn*CHAINLEN;ii+=CHAINLEN) { /* nullize pointers */
	  *(PSTCAST[ii].back) = (SynS *)NULL; /* remove old pre-post pointer */
        }      
        free(PSTCAST);
        free(QPRCAST);  
      }
      nsyn = maxsyn = qlen = 0.;
      lpst = queu = qptr = 0.; 
  }
  ENDVERBATIM
}

PROCEDURE init_arrays(num) {
  VERBATIM {
    int ii, msn;
    if (_lnum == maxsyn) {
      printf("Clearing array\n");
      for (ii=0;ii<nsyn*CHAINLEN;ii+=CHAINLEN) { /* nullize pointers */
	*(PSTCAST[ii].back) = (SynS *)NULL; /* remove old pre-post pointer */
      }
      nsyn = 0;
    } else {
      maxsyn = (double)((int)_lnum);
      msn = maxsyn*CHAINLEN;
      qlen = SPKSTORE*msn;

      if (lpst == 0) {
        nsyn = 0.;
        lpst = (double)((unsigned long)ecalloc(msn,(sizeof(SynS)))); 
        queu = (double)((unsigned long)ecalloc(((int)qlen), (sizeof(QueU)))); 
      } else {
        if (nsyn > maxsyn) {
          printf("Shrinking array\n");
          for (ii=maxsyn;ii<nsyn*CHAINLEN;ii+=CHAINLEN) { /* nullize pointers */
	    /* remove old pre-post pointer */
	    *(PSTCAST[ii].back) = (SynS *)NULL;
          }
          nsyn = maxsyn;
        } else {
          printf("Expanding array %09d\n",QPRCAST->cpost);
        }
        lpst = (double)((unsigned long)erealloc((PSTCAST),((int)(msn*sizeof(SynS))))); 
        queu = (double)((unsigned long)erealloc((QUECAST),((int)(qlen*sizeof(QueU))))); 
        for (ii=0;ii<(int)nsyn*CHAINLEN;ii+=CHAINLEN) {
          *(PSTCAST[ii].back) = &(PSTCAST[ii]); /* restore old pre pointers */
        }
      }

      /* printf("Initializing arrays: Syn= %d bytes,Queue= %d bytes.\n",
	 (int)msn*sizeof(SynS),(int)qlen, sizeof(QueU)); */
      QPRCAST->qq = QUECAST;
    }
  }
 ENDVERBATIM

}

VERBATIM
static void hshake();
ENDVERBATIM
: 2 arguments - index, presyn link
PROCEDURE setlink() {
VERBATIM { 
  int ii, x, new;
  SynS *sns;
  double ptemp;
  PreL *ppsyn;
    
  if (ifarg(2)) { 
    x = (int)*getarg(1);
    ptemp = *getarg(2);		/* pick up pointer */
  } else {
    x = nsyn;
    ptemp = *getarg(1);
  }

  if (x >= maxsyn) {
    init_arrays(maxsyn+POSTINC); /* #DEFINE POSTINC 5 */
  } 
  if (x > nsyn) {
    hoc_execerror("Can't leave empty pointers: see nsyn for current array index.", 0);
  }

  sns = &(PSTCAST[x*(int)CHAINLEN]);	/* postsynaptic entry */
  if ((ppsyn = (PreL *)((unsigned long)ptemp)) == (PreL *)NULL) { /* presyn loc */
    hoc_execerror("Presyn not initialized.", 0); }
  if (ppsyn->link2 != ptemp) {
    hoc_execerror("Invalid (changed) link.", 0); }

  if (x == nsyn) {
    nsyn++;		/* a new entry */
    new = 1;
  } else { 
    /* should generate error if try to change something in middle of chain */
    if (sns->chainlen == -2) {
      hoc_execerror("Internal error: Index used must be multiple of CHAINLEN.", 0);
    }
    new = 0; 
  }

  hshake(sns,ppsyn,new);

  x *= (int)CHAINLEN;
  for (ii=x;ii < x + CHAINLEN;ii++) {
    sns = &(PSTCAST[ii]);	/* access the entry */
    /* initialize postsyn stuff to reasonable values */    
    /* delay, gmax and all codes should be initialized separately */
    sns->del = DELAY;
    sns->ucode = -1;
    sns->chainlen = (ii==x)?CHAINLEN:-2;
    sns->pgm = 1.0;

    /* initialize the presynaptic stuff */
    /* Rcurr and last will be initialized with init */
    sns->index = x;
    sns->qpt = QPRCAST; 
  }
}
ENDVERBATIM
}

: manipulate the presynaptic list remotely
VERBATIM
static void hshake(ss,pl,flag) 
/* ls will be a pointer to presyn cell's array of pointers */
SynS *ss;
PreL *pl;
int flag; /* flag == 1 if this is a brand new entry */
{
  int ii;
  double *nn;
  double *mx;

  nn = pl->nnpre; /* pointer to presyn npre */
  mx = pl->mxpre; /* pointer to presyn maxpre */

  /* erase presyn pointer if this has been set before */
  if (flag == 0) {		/* an old entry */
    /* fall out of loop if a pointer exists already */
    for (ii=0;ii<(*nn) && (pl->plst[ii])!=ss;ii++); /* no body */
    if (ii < (*nn)) {	/* fell out of loop */
      printf("Maintaining pointers.\n");
    } else {
      flag = 1;			/* now we do have to create a new pointer */
      printf("Erasing pointer (C%09d,Pr%09d,Po%09d).\n",
	     ss->ucode,*(ss->pcpre),QPRCAST->cpost);
      *(ss->back) = (SynS *)NULL; /* remove old pre-post pointer */
    }
  }

  if (flag == 1) {		/* create a new pointer */
    /* if necessary, allocate space for the pointer system */
    if ((*nn) == 0.) { 
      (*mx) = PREINC;  /* #DEFINE PREINC 50 */
      pl->plst = (SynS **)ecalloc(((int)(*mx)), sizeof(SynS *)); 
    } else if ((*nn) == (*mx)) { /* need to create more space */
      (*mx) += PREINC;
      pl->plst = (SynS **)realloc(pl->plst,((int)(*mx)) * sizeof(SynS *)); 
      /* reassign all back pointers */
      for (ii=0;ii<(*nn);ii++) {
	(pl->plst[ii])->back = &(pl->plst[ii]);
      }
    }

    /* the handshake */
    pl->plst[(int)(*nn)] = ss;	/* pre  -> post */
    ss->back = &(pl->plst[(int)(*nn)]); /* post -> pre */
    ss->pcpre = &(pl->cpre); /* post -> precode */
    (*nn)++;
  }
}
ENDVERBATIM

: 1 or 2 args, get bzw. set user code
FUNCTION code() {
  VERBATIM { int ii;
    if (ifarg(1)) {
      ii = (int)*getarg(1);
      if (ii < 0) { ii = nsyn+ii; }
      if (ii >= nsyn || ii < 0) {hoc_execerror("array index out of bounds", 0);}
      if (ifarg(2)) {  /* look for a second arg to do a set */
	(PSTCAST[ii*(int)CHAINLEN]).ucode = (int)*getarg(2); }
      _lcode = (PSTCAST[ii*(int)CHAINLEN]).ucode;
    } else {
      for (ii=0;ii<nsyn*CHAINLEN;ii+=CHAINLEN) {
	  printf("%d,%09d  ",ii,(PSTCAST[ii]).ucode);
	}
      _lcode = 1.0;
    }}      
  ENDVERBATIM
}
 
: 2 args => set postsyn code
FUNCTION post() { 
  VERBATIM { 
    if (ifarg(2)) { QPRCAST->cpost = (int)*getarg(2); }
    _lpost = QPRCAST->cpost;
  }
  ENDVERBATIM
}
 
: 1 or 2 args, get bzw. set delay
FUNCTION delay() {
  VERBATIM { int ii,jj;
    if (ifarg(1)) {
      ii = (int)*getarg(1);
      if (ii < 0) { ii = nsyn+ii; }
      if (ii >= nsyn || ii < 0) {hoc_execerror("array index out of bounds", 0);}
      if (ifarg(2)) {  /* look for a second arg to do a set */
	for (jj=ii*CHAINLEN;jj<(ii+1)*CHAINLEN;jj++) {
 	  (PSTCAST[jj]).del = *getarg(2) + DELAY; }}
      _ldelay = (PSTCAST[ii*(int)CHAINLEN]).del;
    } else {
      for (ii=0;ii<nsyn*CHAINLEN;ii+=CHAINLEN) {
	  printf("%d,%g  ",ii,(PSTCAST[ii]).del);
	}
      _ldelay = 1.0;
    }}
  ENDVERBATIM
}

: 1 or 2 args, get bzw. set individual percent gmax
: user expects to see a gmax in uS, but what is stored is a 
: fraction of the total gmax (a global variable)
: therefore can turn off all GABAA with gmax_GABAA = 0
FUNCTION gmax() { 
  VERBATIM { int ii,jj;
    if (ifarg(1)) {
      ii = (int)*getarg(1);
      if (ii < 0) { ii = nsyn+ii; }
      if (ii >= nsyn || ii < 0) {hoc_execerror("array index out of bounds", 0);}
      if (ifarg(2)) {  /* look for a second arg to do a set */
 	for (jj=ii*CHAINLEN;jj<(ii+1)*CHAINLEN;jj++) {
 	  (PSTCAST[jj]).pgm = *getarg(2); }}
      _lgmax = (PSTCAST[ii*(int)CHAINLEN]).pgm;
    } else {
      printf("Multiply by %g to get effective gmax.\n",GMAX);
      for (ii=0;ii<nsyn*CHAINLEN;ii+=CHAINLEN) {
	  printf("%d,%g  ",ii,(PSTCAST[ii]).pgm);
	}
      _lgmax = 1.0;
    }}
  ENDVERBATIM
}

: called with 0,1 or 2 args
: 0 => print out info about all the presyns
: 1 (index) => use as index to return this pre code
:              -1 => print out pointer structure
: 2 (index,code )=> set precode for this index
FUNCTION pre() {
  VERBATIM { int x,ii;
    x = -2;  /* flag -> -1 then print out pointers */
    if (ifarg(1)) { x = (int)*getarg(1); }
    if (x >= 0) {
      if (x >= nsyn) { hoc_execerror("array index out of bounds", 0);}
      x *= CHAINLEN;
      _lpre = *((PSTCAST[x]).pcpre);
    } else {	  
      for (ii=0;ii<nsyn*CHAINLEN;ii+=CHAINLEN) {
	printf("%1d Del:%-2g Gmax: %-4g Cd:%09d Pre:%09d Post:%09d\n",
	       (PSTCAST[ii]).index,
	       (PSTCAST[ii]).del,
	       (PSTCAST[ii]).pgm,
	       (PSTCAST[ii]).ucode,
	       *((PSTCAST[ii]).pcpre),
	       QPRCAST->cpost);
	if (x == -1 || &(PSTCAST[ii]) != *((PSTCAST[ii]).back)) {
	  printf("\t\t%s: %x -> %x\n",
		 ((&(PSTCAST[ii])==*((PSTCAST[ii]).back))?"OK":"POINTER MISMATCH"),
		 &(PSTCAST[ii]),
		 (PSTCAST[ii]).back);
	}
      }
      _lpre = 1.0;
    }
  } 
  ENDVERBATIM
}

: return link for matching with presyn
FUNCTION link(ii) { 
  VERBATIM { 
    if (_lii >= nsyn) { hoc_execerror("array index out of bounds", 0);}
    _llink = (double)(unsigned long)(PSTCAST[(int)_lii]).pcpre;
  }
  ENDVERBATIM
}

FUNCTION check() {
  VERBATIM { 
    int ii;
    for (ii=0; ii<nsyn*CHAINLEN && _lcheck==1.; ii+=CHAINLEN) {
      if (&(PSTCAST[ii]) != *((PSTCAST[ii]).back)) {
	printf("****************************************************************\n");
	printf("ERROR:: Index:%3d,C%09d,Pr%09d,Po%09d, Delay:%6g (%x->%x->%x)\n",
	       (PSTCAST[ii]).index,
	       (PSTCAST[ii]).ucode,
	       *((PSTCAST[ii]).pcpre),
	       QPRCAST->cpost,
	       (PSTCAST[ii]).del,
	       &(PSTCAST[ii]),
	       (PSTCAST[ii]).back,
	       *((PSTCAST[ii]).back));
	printf("****************************************************************\n");
	_lcheck = -1.;
      }
    }
    _lcheck = nsyn;  /* all pointers have to be active */
  }
  ENDVERBATIM
}

PROCEDURE initq() { 
  VERBATIM {
  int ii;
  QPRCAST->dead = Cdur + Deadtime;
  begsyn = endsyn = newspk = 0.;
  for (ii=0;ii<qlen;ii++) {
    QUECAST[ii].time = 1.e20;
    QUECAST[ii].index = -1;
  }
  QPRCAST->qterm = -1.e2;
  for (ii=0;ii<nsyn*CHAINLEN;ii++) {
    PSTCAST[ii].Rcurr = 0.;
    PSTCAST[ii].last = -1.e2;
    PSTCAST[ii].spkt = -1.e2;
    if (PSTCAST[ii].chainlen > 1) { PSTCAST[ii].chainptr = -1; }
  }
  }
  ENDVERBATIM
}

: begsyn is the top of the line when the spike is first utilized to generate synapse
: updates the value in the queue so endsyn can detect the end of the syn pulse
: doesn't return value (since in a struct)
PROCEDURE popqh1(aug) { 
  VERBATIM {  
    if (QUECAST[(int)begsyn].time == 1e20) {
      printf("%x %g %g ",QUECAST,newspk,qlen);
      hoc_execerror("Error: queue exhausted.\n",0);
    } else {  /* augment the time by Cdur */
      QUECAST[(int)begsyn].time += _laug;
      begsyn++;
      if (begsyn == qlen) { begsyn = 0.; }
    }
  }
  ENDVERBATIM
}
 
: endsyn is queried for the termination time of the synaptic square wave
PROCEDURE popqh2() { 
  VERBATIM {  
    QUECAST[(int)endsyn].time = 1.e20;  /* clear the entry */
    endsyn++;
    if (endsyn == qlen) { endsyn = 0.; }
  }
  ENDVERBATIM
}

: DEBUGGING ROUTINES

FUNCTION getdbx(c,x) {
VERBATIM { 
  switch ((int)_lc) {
  case 1:
    _lgetdbx = (PSTCAST[(int)_lx]).last;
    break;
  case 2:
    _lgetdbx = (PSTCAST[(int)_lx]).Rcurr;
    break;
  case 3:
    _lgetdbx = (PSTCAST[(int)_lx]).spkt;
    break;
  case 4:
    _lgetdbx = (QUECAST[(int)_lx]).time;
    break;
  case 5:
    _lgetdbx = (QUECAST[(int)_lx]).index;
    break;
  case 6:
    _lgetdbx = QPRCAST->qterm;
    break;
  case 7:
    _lgetdbx = QPRCAST->dead;
    break;
  case 8:
    _lgetdbx = (PSTCAST[(int)_lx]).chainlen;
    break;
  case 9:
    _lgetdbx = (PSTCAST[(int)_lx]).chainptr;
    break;
  default:
    hoc_execerror("UNAVAILABLE IN DBX",0);
    _lgetdbx = -1;
    break;
  }}
  ENDVERBATIM
}

PROCEDURE prq() {
  VERBATIM { int ii;
    printf("new:%g beg:%g end:%g /%g\n",newspk,begsyn,endsyn,qlen);
    for (ii=endsyn;ii!=newspk;ii=((ii==qlen-1)?0:ii+1)) {
        printf("%d %8g%8d\n",ii,(QUECAST[ii]).time,(QUECAST[ii]).index);
    }
  }
  ENDVERBATIM
}


