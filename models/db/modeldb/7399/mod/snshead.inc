: $Id: snshead.inc,v 1.22 1998/01/07 20:22:30 billl Exp $

VERBATIM

#define PSTCAST ((SynS *)((unsigned long) lpst))
#define QUECAST ((QueU *)((unsigned long) queu))
#define PRECAST ((PreL *)((unsigned long) link))
#define QPRCAST ((QptR *)((unsigned long) qptr))

#define PREINC  50.  /* size increment for growing the presyn array */
#define POSTINC 5.   /* size increment for growing the postsyn array */

typedef struct QueU { /* postsynaptic queue */
  double time;  /* records spike time + delay */
  short index;  /* records index of associated SynS struct */
} QueU;

/* holds all the info for accessing the queu and its tail */
typedef struct QptR {
  struct QueU *qq;     /* will point to the queu */
  double *nspk; /* will point to newspk */
  double *head; /* will point to begsyn */
  double *qln;  /* points to qlen */
  double dead;  /* Cdur + Deadtime */
  double qterm;  /* the last time put on this queue */
  int cpost;  /* user code for info about postsyn cell */
} QptR;

typedef struct SynS { /* postsynaptic structure */
  /* this stuff used presynaptically */
  short index;  /* index of this item */
  double del;   /* param: individualized delay set by user */
  double spkt;  /* state: time when next spike can begin */
  struct QptR *qpt;   /* indirection to get at the queu */

  /* for allowing mult spikes to add up */
  short chainlen;
  short chainptr;

  /* user codes */
  int ucode;  /* user code for user convenience, an arbitrary index */
  /* this stuff used postsynaptically */
  double pgm;   /* percent of gmax - set by user */
  double Rcurr; /* holds current value of R (release status), internal */
  double last;  /* time last syn will end, internal */
  int *pcpre;   /* pointer to access the presyn code */
  struct SynS **back;  /* pointer to presyn pointer to this entry */
} SynS;

/* holds all the info for accessing presyn list */
typedef struct PreL {
  int cpre;   /* user code for info about presyn cell */
  struct SynS **plst;        /* will point to post for postsyn nrn*/
  double *nnpre;     /* number of pointers currently there */
  double *mxpre;     /* number permitted before realloc required */
  double link2;       /* a copy of the link information to confirm correct hook up */
} PreL;

ENDVERBATIM
