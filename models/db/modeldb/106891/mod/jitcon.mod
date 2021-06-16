: $Id: jitcon.mod,v 1.20 2007/12/28 02:58:37 billl Exp $

:* main COMMENT
COMMENT
simulation to accompany "Just in time connectivity for large spiking networks"
jitcon.mod code extracted from intf.mod627 hence retains a small amount of nonfunctional code
  some of this has been purposely left in there since might be of use in future versions with
  greater functionality
current version performs JitEvent with callbacks() and JitWeights but does not perform JitCon --
  ie does not recreate the divergence vector; instead divergence vector for all cells is created
  with a calls from each JitCon to setdvi()
JitCon is a presynaptic process which should be connected by a NetCon to a compartmental soma eg
  ncn=new NetCon(jitcon, precell.newsyn(type,ident)) // eg bsticknet.hoc_35:174
in order to simplify establishing connectivity, the postsynaptic expsyn (from myexpsyn.mod)
  maintains information about cell ID and cell type of its cell
because of the difficulty of following pointers through hoc-level templates, all required
  information is maintained in parallel lists which can be readily accessed in nmodl:
    cells is a list of cells (not used by JitCon)
    ce is a list of associated JitCon's
    aml is list of lists of AMPA synapses
    gal is list of lists of GABAA synapses
access to all of these lists as well as to numerous shared data structures (eg weight array,
  connectivity arrays (pmat and div), delay array, etc. is achieved through 
  calls to the hoc level performed in jitcondiv(): eg HPTR("wmat") returns a pointer to data
  structure 'wmat' in hoc, HVAL("allcells") returns the value of param 'allcells' in hoc
JitCon function is controlled by the 'jcn' flag which can be set or cleared in hoc
  jcn==0 -- JitCon NetReceive block calls net_event() which will activated associated NetCons
            ie JitCon is not being used
  jcn==1 -- calls jitcon() which places the first callback on the queue
the callback is a self_event() which uses descending sequential flag<0 to indicate position in
  the divergence vector; NET_RECEIVE() with flag<0 calls callback()
  callback() places another callback() self_event on the queue and also makes a direct call
  to the NET_RECEIVE block of the indicated expsyn
Additional versions of JitCon were available in intf.mod by using jcn=2 or jcn=3; these have
  not been re-implemented here though the code is included for possible future use.
  In general the difficulty in all of these implementations lies in precisely recreating the
  random values obtained via hoc calls -- identical NetCon/JitCon networks are useful for network
  validation
ENDCOMMENT

:* main VERBATIM block
VERBATIM
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h> /* contains LONG_MAX */

extern void* vector_arg();
extern double hoc_call_func(Symbol*, int narg);
extern double* hoc_pgetarg();
extern FILE* hoc_obj_file_arg(int narg);
extern int list_vector_px(Object *ob, int i, double** px);
extern int list_vector_px2 (Object *ob, int i, double** px, void** vv);
extern Object** hoc_objgetarg();
extern int ivoc_list_count(Object*);
extern Object* ivoc_list_item(Object*, int);
extern Symbol *hoc_get_symbol();
extern Symbol *hoc_lookup();
extern Point_process* ob2pntproc(Object*);
extern double mcell_ran4();
extern int hoc_is_double_arg(int narg);
static void hxe() { hoc_execerror("",0); }
#if defined(t)
static void initmodel();
#else
static initmodel();
#endif
extern int stoprun;
extern double hoc_epsilon;
extern short *nrn_artcell_qindex_;
extern double nrn_event_queue_stats();
extern void clear_event_queue();
extern double *vector_newsize();
extern Objectdata *hoc_objectdata;
extern int nrn_mlh_gsort();
extern int cmpdfn();
extern unsigned int *scrset();

#define CTYN 2  // number of cell types being used
#define PI 3.141592653589793115997963468544
#define nil 0
#define SOP (((id0*) _p_sop)->vp)
#define IDP (*((id0**) &(_p_sop)))
#define NSW 20  // just store voltages
#define NSV 7  // 6 state variables (+ 1 for time)
#define FOFFSET 100 // flag offset for NET_RECEIVE()
#define DELMIN 1e-5 // 10 ns is minimum delay to bother using queue -- otherwise consider sync
#define DELM(X,Y) (*(pg->delm+(X)*CTYPi+(Y)))
#define DELD(X,Y) (*(pg->deld+(X)*CTYPi+(Y)))
#define DVG(X,Y) ((int)*(pg->dvg+(X)*CTYPi+(Y)))
#define WMAT(X,Y) (*(pg->wmat+(X)*CTYPi+(Y)))
#if defined(t)
#define HVAL(X) (*(hoc_objectdata[(hoc_get_symbol((X)))->u.oboff]._pval))
#define HPTR(X) (hoc_objectdata[(hoc_get_symbol((X)))->u.oboff]._pval)
#else
#define HVAL(X) (*(hoc_objectdata[(hoc_get_symbol((X)))->u.oboff].pval))
#define HPTR(X) (hoc_objectdata[(hoc_get_symbol((X)))->u.oboff].pval)
#endif

typedef struct POSTGRP { // postsynaptic group
  double *dvg; double *delm; double *deld; double *ix; double *ixe; double *wmat;
  struct POSTGRP *next;
} postgrp;

typedef struct ID0 {
  postgrp *pg;
  Point_process **dvi; // each cell has a divergence list
  double *del;        // each syn has its own intrinsic delay
  unsigned char *sprob;    // each syn has a firing probability 0-255->0-1
  unsigned int dvt;
  unsigned int  id;
  // type -> jcn MUST REMAIN unbroked BLOCK -- see flag()
  // when adding flags also augment iflags, iflnum
  // only use first 3 letters with flag() -- see iflags
  unsigned char     type;  // | 
  unsigned char     col;   // |
           char     dbx;   // |
  unsigned char     jcn;   // |
  // end BLOCK
} id0;

// globals -- range vars must be malloc'ed in the CONSTRUCTOR
static id0 *ip, *qp;
static postgrp *pg;
static Object *ce, *gal, *aml;
static Point_process *pmt, *tpnt;
static char *name;
// iflags string use to find flags -- note that only 1st 3 chars are used to identify
static char iflags[100]="typ col dbx jcn"; 
static char iflnum=4;      // turn on after generating an error message
static double *lop(); // accessed by all JitCon
static unsigned int cesz; // number of cells
static unsigned int sead; // 'sead' vs global 'seed'
static int  CTYPi, STYPi, scrsz; // from labels.hoc
static double allcells, qlimit, *scr;
void*    ww[NSW];
static int AM=0, NM=1, GA=2, GB=3, SU=3, IN=4, DP=2; // from labels.hoc
static char* CNAME[10];
// static int cty[CTYN]={3, 4}; // SU IN -- can only use constants here
static int cty[CTYN], process;
static double wts[1];  // for callback to use as the wts[] pointer
ENDVERBATIM

:* NEURON, PARAMETER, ASSIGNED blocks
NEURON {
  ARTIFICIAL_CELL JitCon
  POINTER sop                          :::: Structure pointer for C-level range vars
}

PARAMETER {
  sop=0
}

ASSIGNED {
}

:* CONSTRUCTOR, DESTRUCTOR, INITIAL
:** CONSTRUCT: create a structure to save the identity of this unit and char integer flags
CONSTRUCTOR {
  VERBATIM 
  { int lid,lty,lco,i;
    if (ifarg(1)) { lid=(int) *getarg(1); } else { lid= UINT_MAX; }
    if (ifarg(2)) { lty=(int) *getarg(2); } else { lty= -1; }
    if (ifarg(3)) { lco=(int) *getarg(3); } else { lco= -1; }
    _p_sop = (void*)ecalloc(1, sizeof(id0)); // important that calloc sets all flags etc to 0
    ip = IDP;
    ip->id=lid; ip->type=lty; ip->col=lco; ip->pg=0x0; ip->dvi=0x0; ip->sprob=0x0;
    ip->jcn = 0;
    process=(int)getpid();
    CNAME[SU]="SU"; CNAME[DP]="DP"; CNAME[IN]="IN";
  }
  ENDVERBATIM
}

DESTRUCTOR {
  VERBATIM { 
  free(IDP);
  }
  ENDVERBATIM
}

:** INITIAL
INITIAL {
}

:* NET_RECEIVE
NET_RECEIVE (w) { LOCAL tmp,jcn
  VERBATIM
  ip = IDP;
  _ljcn=ip->jcn;
  tpnt = _pnt; // this pnt
  // jitcon off -- use NetCon
  // the C version of mod net_event call; pass on via regular NetCon connections 
  if (! ip->jcn) { 
    net_event(tpnt, t); 
  // else jitcon on -- if flag==0 (external event) -- initiate sequence
  } else if (_lflag==0) { 
    jitcon(t);
  // jitcon on -- flag<0 - a callback
  } else if (_lflag<0) { 
    callback(_lflag); 
  }
  ENDVERBATIM
}

:* ancillary functions
:** jitcon() creates divergence and delays from rand seed
: jcn flags:
: 0 NetCons                                   jcn==0
: 1 Jitcon with callback with stored pointers jcn==1
: 2 Jitcon with Jitcon (regenerated pointers) jcn==2 // not re-implemented
: 3 Jitcon without JitEvent                   jcn==3 // not re-implemented
PROCEDURE jitcon (tm) {
  VERBATIM {
  double mindel, randel, idty, *x; int prty, poty, i, j, k, dv, dvt; 
  Point_process *pnt; void* voi;
  // qsz = nrn_event_queue_stats(stt);
  // if (qsz>=qlimit) { printf("qlimit %g exceeded at t=%g\n",qlimit,t); qlimit*=2; }
  ip=IDP;
  prty=(int)ip->type;
  if (ip->jcn==1) { //
    if (!pg) {printf("No network defined -- must run jitcondiv()\n"); hxe();}
#if defined(t)
    if (ip->dvt>0) net_send((void**)0x0, wts,tpnt,ip->del[0]+t,-1.); // first callback
#else
    if (ip->dvt>0) net_send((void**)0x0, wts,tpnt,ip->del[0],-1.); // first callback
#endif
  } else if (ip->jcn==3) { mkdvi(); // make divergence lists on the fly
  } else if (ip->jcn==2) { // true JitCon: not currently working
    // concept here is to replicate the code in bsticknet.hoc_35:207 -- connec()
    // minor difficulty is in recreating rdmuniq() in C
    // a dirty, trick would be to just call the hoc routine to set this
    // I have not bothered to do this since it does not seem like something anyone
    // would actually want to use
    for (i=0,k=0,dvt=0;i<CTYN;i++) { // calculated total divergence for malloc'ing
      poty=cty[i];
      printf("%d:%d:%d ",prty,poty,(int)DVG(prty,poty));
      dvt+=(int)DVG(prty,poty);
    }
    for (i=0,k=0,dvt=0;i<CTYN;i++) { // cell types in cty[]
      if (dv>0) {
        mcell_ran4(&sead, &randel , 1, 2.*DELD(prty,poty)); // delays
        randel+=((_ltm-t)+DELM(prty,poty)-DELD(prty,poty));
        if (randel<=0) randel= -randel;
        if (dv>scrsz) {
          printf("A:Divergence exceeds scrsz: %d>%d for %d->%d\n",dv,scrsz,prty,poty); hxe(); }
        mcell_ran4(&sead, scr ,  dv, pg->ixe[poty]-pg->ix[poty]+1); // indices
        if (ifarg(2)) { // dump -- not needed since can replicate the sequence from hoc
          voi=vector_arg(2);
          x=vector_newsize(voi,k+dv+1);
          x[k++]= -randel; // store as a negative number
          for (j=0;j<dv;j++,k++) x[k]=floor(scr[j]+pg->ix[poty]);
        } else for (j=0;j<dv;j++) {
          pnt=(Point_process *)(ivoc_list_item(ce, (int)(scr[j]+pg->ix[poty])))->u.this_pointer;
          idty=(double)(FOFFSET+ip->id)+0.1*(double)ip->type+0.01;
#if defined(t)
          net_send((void**)0x0, wts, pnt, randel+t, idty);
#else
          net_send((void**)0x0, wts, pnt, randel, idty);
#endif
        }
      }
    }
  } 
  }
  ENDVERBATIM  
}

PROCEDURE callback (fl) {
  VERBATIM {
  int ii,jj; double idty, del; double weed, prid, prty, poid, poty, w;
  unsigned int valseed;
  ii=(unsigned int)((-_lfl)-1); // -1,-2,-3 -> 0,1,2
  ip=IDP;
  idty=(double)(FOFFSET+ip->id)+0.1*(double)ip->type+0.01;
  jj=ii+1;
  if (jj<ip->dvt) {
    del= ip->del[jj] - ip->del[ii];
#if defined(t)
    net_send((void**)0x0, wts,tpnt,del+t,(double) -(jj+1)); // next callback
#else
    net_send((void**)0x0, wts,tpnt,del,(double) -(jj+1)); // next callback
#endif
  }
  if (ip->sprob[ii]) { // keep this for now in order to handle first callback 
    // 'idty' was meant to be sent as a flag so postsyn would have preid to generate its weights
    // this isn't necessary -- can generate the weights right here
    poid=ip->dvi[ii]->_prop->param[3]; // BAD -- absolute addresses used here
    poty=ip->dvi[ii]->_prop->param[4];
    prid=ip->dvi[ii]->_prop->param[5];
    prty=ip->dvi[ii]->_prop->param[6];
    if ((((double)ip->type)!=prty) || (((double)ip->id)!=prid)) {
      printf("callback() ERR: %g %g %g %g %g %g\n",\
             ((double)ip->type),prty,((double)ip->id),prid,poty,poid); hxe(); }
    weed=prty*allcells+poty*100+prid*10+poid;
    // vw1.setrnd(4,2*0.01*w,weed) vw.add(0.99*w) // from (bsticknet.hoc_32:235)
    valseed=(unsigned int)weed; w=WMAT((int)prty,(int)poty);
    mcell_ran4(&valseed, &wts, 1, 2*0.01*w); // generate 1 value
    // printf("%g %g %g %g %g\n",w,wts[0],weed,prid,poid);
    wts[0]+=0.99*w; // note that weight is created presynaptically rather than postsynaptically
                    // as in intf.mod
    (*pnt_receive[ip->dvi[ii]->_prop->_type])(ip->dvi[ii], wts, 0.0); } 
  } 
  ENDVERBATIM
}

:* mkdvi() -- create random graph on the fly -- not updated for jitcon.mod
PROCEDURE mkdvi () {
VERBATIM {
  int i,j,k,prty,poty,dv,dvt,dvii; double *x, *db, *dbs; 
  Object *lb;  Point_process *pnnt, **da, **das;
  ip=IDP; ip->pg=pg; // this should be called right after jitcondiv()
  prty=ip->type;
  sead=((unsigned int)ip->id)*1e6;
  for (i=0,k=0,dvt=0;i<CTYN;i++) { // dvt gives total divergence
    poty=cty[i];
    dvt+=DVG(prty,poty);
  }
  da =(Point_process **)malloc(dvt*sizeof(Point_process *));
  das=(Point_process **)malloc(dvt*sizeof(Point_process *)); // das,dbs for after sort
  db =(double *)malloc(dvt*sizeof(double)); // delays
  dbs=(double *)malloc(dvt*sizeof(double)); // delays
  for (i=0,k=0,dvii=0;i<CTYN;i++) { // cell types in cty[]
    poty=cty[i];
    dv=DVG(prty,poty);
    if (dv>0) {
      sead+=dv;
      if (dv>scrsz) {
        printf("B:Divergence exceeds scrsz: %d>%d for %d->%d\n",dv,scrsz,prty,poty); hxe(); }
      mcell_ran4(&sead, scr ,  dv, pg->ixe[poty]-pg->ix[poty]+1);
      for (j=0;j<dv;j++) {
        if (!(lb=ivoc_list_item(ce,(unsigned int)floor(scr[j]+pg->ix[poty])))) {
          printf("JitCon:callback %g exceeds %d for list ce\n",floor(scr[j]+pg->ix[poty]),cesz); 
          hxe(); }
        pnnt=(Point_process *)lb->u.this_pointer;
        da[j+dvii]=pnnt;
      }
      mcell_ran4(&sead, scr , dv, 2*DELD(prty,poty));
      for (j=0;j<dv;j++) {
        db[j+dvii]=scr[j]+DELM(prty,poty)-DELD(prty,poty); // +/- DELD
        if (db[j+dvii]<0) db[j+dvii]=-db[j+dvii];
      }
      dvii+=dv;
    }
  }
  gsort2(db,da,dvt,dbs,das);
  ip->del=dbs;   ip->dvi=das;   ip->dvt=dvt;
  ip->sprob=(unsigned char *)malloc(dvt*sizeof(char *)); // release probability
  for (i=0;i<dvt;i++) ip->sprob[i]=1; // start out with all firing
  free(da); free(db); // keep das,dbs which are assigned to ip->dvi bzw ip->del
  }
ENDVERBATIM
}

:* getdvi() retrieves divergence and delay information after setdvi()
FUNCTION getdvi () {
  VERBATIM 
  {
  int j,dvt; double *dbs, *x;
  void* voi; Point_process **das;
  ip=IDP; ip->pg=pg; // this should be called right after jitcondiv()
  dvt=ip->dvt;
  dbs=ip->del;   das=ip->dvi;
  _lgetdvi=(double)dvt; 
  if (!ifarg(1)) return _lgetdvi; // just return the divergence value
  if (hoc_is_double_arg(1)) { // return the effective divergence value
    
  }
  voi=vector_arg(1);
  x=vector_newsize(voi,dvt);
  for (j=0;j<dvt;j++) {
    x[j]=(double)das[j]->_prop->param[3]; // BAD -- no way to check that this is in fact "id"
  }
  voi=vector_arg(2);
  x=vector_newsize(voi,dvt);
  for (j=0;j<dvt;j++) x[j]=dbs[j];
  if (ifarg(3)) {
    voi=vector_arg(3);
    x=vector_newsize(voi,dvt);
    for (j=0;j<dvt;j++) x[j]=(double)ip->sprob[j];
  }
  }
ENDVERBATIM
}

:* setdvi() is called with postyn id vector and delay vector to set up connections
PROCEDURE setdvi () {
VERBATIM {
  int i,j,k,dvt,dvu,ddvi,lbcnt; double *y, *db, *dbs, id;
  void* voi; Object *lb, *lc, *syo; Point_process *pnnt, **da, **das;
  ip=IDP; ip->pg=pg; id=(double)ip->id;
  dvt=vector_arg_px(1, &y); // dvt is divergence
  i=vector_arg_px(2, &db);
  if (i != dvt) {printf("setdvi() ERR vec sizes: %d %d\n",dvt,j); hxe();}
  if (ip->type==IN) syo=gal; else syo=aml; // list depends on presyn type
  da =(Point_process **)malloc(dvt*sizeof(Point_process *));
  das=(Point_process **)malloc(dvt*sizeof(Point_process *)); // das,dbs for after sort
  dbs=(double *)malloc(dvt*sizeof(double)); // delays
  for (j=0,i=0;j<dvt;j++) {
    lb=ivoc_list_item(syo,(unsigned int)y[j]); // lb is another list
    if (!lb) { printf("JitCon:callback %g exceeds %d for list ce\n",y[j],cesz); hxe(); }
    lbcnt=ivoc_list_count(lb);
    for (k=0;k<lbcnt;k++) {
      lc=ivoc_list_item(lb,k);
      if (((Point_process *)lc->u.this_pointer)->_prop->param[5] == id) {
        da[j]=(Point_process *)lc->u.this_pointer;
        break;
      }
    }
    if (k==lbcnt) { printf("setdvi() ERR: id %g not found\n",id); hxe(); }
  }
  gsort2(db,da,dvt,dbs,das);
  ip->del=dbs;   ip->dvi=das;   ip->dvt=dvt;
  ip->sprob=(unsigned char *)malloc(dvt*sizeof(char *)); // release probability
  for (j=0;j<dvt;j++) ip->sprob[j]=1; // start out with all firing
  free(da);
  }
ENDVERBATIM
}

: prune(p[,rand_seed]) // prune synapses with prob p [0,1], ie 0.1 prunes 10% of the divergence
: prune(vec) // fill in the pruning vec with binary values from vec
PROCEDURE prune () {
  VERBATIM 
  {
  double *x, p; int nx,j;
  ip=IDP; ip->pg=pg;
  if (hoc_is_double_arg(1)) {
    p=*getarg(1);
    if (p<0 || p>1) {printf("JitCon:pruneERR0:need # [0,1] to prune [ALL,NONE]: %g\n",p); hxe();}
    if (p==1.) printf("JitConpruneWARNING: pruning 100% of cell %d\n",ip->id);
    if (ip->dvt>scrsz) {
      printf("JitConpruneB:Div exceeds scrsz: %d>%d\n",ip->dvt,scrsz); hxe(); }
    for (j=0;j<ip->dvt;j++) ip->sprob[j]=1; // unprune completely
    if (p==0.) return; // now that unpruning is done, can return
    sead=(ifarg(2))?(unsigned int)*getarg(2):(unsigned int)ip->id*1e6;
    mcell_ran4(&sead, scr , ip->dvt, 1.0); // random var (0,1)
    for (j=0;j<ip->dvt;j++) if (scr[j]<p) ip->sprob[j]=0; // prune with prob p
  } else {
    nx=vector_arg_px(1,&x);
    if (nx!=ip->dvt) {printf("JitCon:pruneERRA:Wrong size vector:%d!=%d\n",nx,ip->dvt); hxe();}
    for (j=0;j<ip->dvt;j++) ip->sprob[j]=(unsigned char)x[j];
  }
  }
ENDVERBATIM
}

VERBATIM 
// gsort2() sorts 2 parallel vectors -- delays and Point_process pointers
int gsort2 (double *db, Point_process **da,int dvt,double *dbs, Point_process **das) {
  unsigned int *scr, i;
  scr=scrset(dvt);
  for (i=0;i<dvt;i++) scr[i]=i;
  nrn_mlh_gsort(db, scr, dvt, cmpdfn);
  for (i=0;i<dvt;i++) {
    dbs[i]=db[scr[i]];
    das[i]=da[scr[i]];
  }
}
ENDVERBATIM

PROCEDURE freedvi () {
  VERBATIM
  { 
    int i, poty; id0 *jp;
    jp=IDP;
    if (jp->dvi) {
      free(jp->dvi);
      free(jp->del);
      jp->dvi=0x0;
      jp->del=0x0;
    }
  }
  ENDVERBATIM
}

FUNCTION qstats () {
  VERBATIM {
    double stt[3]; int lct,flag; FILE* tf;
    if (ifarg(1)) {tf=hoc_obj_file_arg(1); flag=1;} else flag=0;
    lct=cty[IDP->type];
    _lqstats = nrn_event_queue_stats(&stt);
    printf("QUEUE: Inserted %g; removed %g\n",stt[0],stt[2]);
    if (flag) {
      fprintf(tf,"QUEUE: Inserted %g; removed %g remaining: %g\n",stt[0],stt[2],_lqstats);
    }
  }
  ENDVERBATIM
}

FUNCTION qsz () {
  VERBATIM {
    double stt[3];
    _lqsz = nrn_event_queue_stats(&stt);
  }
  ENDVERBATIM
}

PROCEDURE qclr () {
  VERBATIM {
    clear_event_queue();
  }
  ENDVERBATIM
}

: intf.jitcondiv() assigns pointers for hoc symbol storage
PROCEDURE jitcondiv () {
  VERBATIM {
  Symbol *sym; int i,j; char name[100];
  ip->pg=(postgrp *)malloc(sizeof(postgrp));
  pg=ip->pg;
  sym = hoc_lookup("ce"); ce = (*(hoc_objectdata[sym->u.oboff].pobj));
  sym = hoc_lookup("aml"); aml = (*(hoc_objectdata[sym->u.oboff].pobj));
  sym = hoc_lookup("gal"); gal = (*(hoc_objectdata[sym->u.oboff].pobj));
  cesz = ivoc_list_count(ce);
  if (cesz!=(i=ivoc_list_count(aml)) || cesz!=(j=ivoc_list_count(gal))) {
    printf("All 3 lists should be same size: ce,aml,gal %d,%d,%d\n",cesz,i,j); hxe(); }
  cty[0]=SU; cty[1]=IN; // set the cell types
  CTYPi=HVAL("CTYPi"); STYPi=HVAL("STYPi"); scrsz=HVAL("scrsz"); allcells=HVAL("allcells");
  pg->ix =HPTR("ix"); pg->ixe=HPTR("ixe"); 
  pg->dvg=HPTR("div"); 
  pg->wmat=HPTR("wmat");
  pg->delm=HPTR("delm"); pg->deld=HPTR("deld");
  scr=HPTR("scr");
  if (!ce) {printf("JitCon jitcondiv ERRA: ce not found\n"); hxe();}
  // make sure no seg error:
  printf("Checking for possible seg error in double arrays: CTYPi==%d: ",CTYPi);
  // can access arbitrary member dvg[a][b] using (&dvg[a*CTYPi])[b] or dvg+a*CTYPi+b
  printf("%d %d %d ",DVG(CTYPi-1,CTYPi-1),(int)pg->ix[CTYPi-1],(int)pg->ixe[CTYPi-1]);
  printf("%g ",WMAT(CTYPi-1,CTYPi-1));
  printf("%g %g ",DELM(CTYPi-1,CTYPi-1),DELD(CTYPi-1,CTYPi-1));
  printf("%d %g\n",scrsz,scr[scrsz-1]); // scratch area for doubles
  }
  ENDVERBATIM  
}

:** probejcd()
PROCEDURE probejcd () {
  VERBATIM {  int i,a[4];
    for (i=1;i<=3;i++) a[i]=(int)*getarg(i);
    printf("CTYPi: %d, STYPi: %d, ",CTYPi,STYPi);
    // printf("div: %d, ix: %d, ixe: %d, ",DVG(a[1],a[2]),(int)ix[a[1]],(int)ixe[a[1]]);
    printf("wmat: %g\n",WMAT(a[1],a[2]));
  }
  ENDVERBATIM  
}

:** vers gives version
PROCEDURE vers () {
  printf("$Id: jitcon.mod,v 1.20 2007/12/28 02:58:37 billl Exp $\n")
}

VERBATIM
//* internal routines
//** lop(LIST,ITEM#) lop and lop set different global variables
// modeled on vector_arg_px(): picks up obj from list and resolves pointers
static double* lop (Object *ob, unsigned int i) {
  Object *lb;
  lb = ivoc_list_item(ob, i);
  if (! lb) { printf("JitCon:lop %d exceeds %d for list ce\n",i,cesz); hxe();}
  pmt=ob2pntproc(lb);
  qp=*((id0**) &((pmt->_prop->dparam)[2])); // #define sop *_ppvar[2].pval
  return pmt->_prop->param;
}

// use stoppo() as a convenient conditional breakpoint in gdb (gdb watching is too slow)
int stoppo () {
}
ENDVERBATIM


: lof can find object information
PROCEDURE lof () {
VERBATIM {
  Object *ob; int num,i,ii,j,k,si,nx;  double *vvo[7], *par; void *vv[7];
  ob = *(hoc_objgetarg(1));
  si=(int)*getarg(2);
  num = ivoc_list_count(ob);
  if (num!=7) { printf("JitCon lof ERR %d>7\n",num); hxe(); }
  for (i=0;i<num;i++) { 
    j = list_vector_px3(ob, i, &vvo[i], &vv[i]);
    if (i==0) nx=j;
    if (j!=nx) { printf("JitCon lof ERR %d %d\n",j,nx); hxe(); }
  }
 }
ENDVERBATIM
}

: ldv prints out divergence arrays
PROCEDURE ldv () {
VERBATIM {
  Object *ob; 
  int dv,num,i,j,ii,prty,poty,nx,offset;  
  double *vvo[100000], *x; void *vv[100000];
  ob = *(hoc_objgetarg(1));
  nx = vector_arg_px(2, &x); // cell# vec
  prty=(int)*getarg(3);
  poty=(int)*getarg(4);
  if (ifarg(5)) offset=(int)*getarg(5); else offset=0;
  num = ivoc_list_count(ob);
  if (num!=nx) {printf("JitCon ldv ERRD %d %d\n",num,nx); hxe(); }
  // if (ix[prty]+offset+num>ixe[prty]) {
  //   printf("JitCon ldv ERR0 %d %d\n",ix[prty]+offset+num,ixe[prty]); hxe(); }
  if (num>1e5) { printf("JitCon ldv ERRA %d>1e5\n",num); hxe(); }
  i=0; nx=list_vector_px3(ob, i, &vvo[i], &vv[i]);
  dv=DVG(prty,poty);
  if (nx!=dv) { printf("JitCon ldv ERRB %d %d\n",dv,nx); hxe(); }
  for (i=1;i<num;i++) { 
    j = list_vector_px3(ob, i, &vvo[i], &vv[i]);    
    if (j!=nx) { printf("JitCon ldv ERRC %d %d\n",j,nx); hxe(); }
  }
  if (ii!=num) printf("INF ldv WARNING: only filled %d of %d columns\n",ii,num);
 }
ENDVERBATIM
}

:** chk([flag]) with flag=1 prints out info on the record structure
:                    flag=2 prints out info on the global vectors
PROCEDURE chk (f) {
  VERBATIM 
  {int i,lfg;
  lfg=(int)_lf;
  ip=IDP;
  if (lfg==1) {
  }
  if (lfg==2) { 
  }
  if (lfg==3) { 
  }
  if (lfg==4) { 
  }
  if (lfg==5) { 
  }}
  ENDVERBATIM
}

:** id() and pid() identify the cell -- printf and function return
FUNCTION pid () {
  VERBATIM 
  printf("JitCon%d(%d/%d@%g) ",IDP->id,IDP->type,IDP->col,t);
  _lpid = (double)IDP->id;
  ENDVERBATIM
}

FUNCTION id () {
  VERBATIM
  if (ifarg(1)) IDP->id = (unsigned int) *getarg(1);
  _lid = (double)IDP->id;
  ENDVERBATIM
}

FUNCTION type () {
  VERBATIM
  if (ifarg(1)) IDP->type = (unsigned char) *getarg(1);
  _ltype = (double)IDP->type;
  ENDVERBATIM
}

FUNCTION col () {
  VERBATIM 
  ip = IDP; 
  if (ifarg(1)) ip->col = (unsigned char) *getarg(1);
  _lcol = (double)ip->col;
  ENDVERBATIM
}

FUNCTION dbx () {
  VERBATIM 
  ip = IDP; 
  if (ifarg(1)) ip->dbx = (unsigned char) *getarg(1);
  _ldbx = (double)ip->dbx;
  ENDVERBATIM
}

:** setting and getting flags: fflag, record,input,jttr
FUNCTION fflag () { fflag=1 }

: flag(name,[val]) set or get the a flag
: seek names from iflags[] and look at location &ip->type -- beginning of flags
: opt 3rd arg to set flag for all of them
FUNCTION flag () {
  VERBATIM 
  {char *sf; int ii,i; unsigned char val;
  ip = IDP;
  sf = gargstr(1);
  for (ii=0;ii<iflnum && strncmp(sf, &iflags[ii*4], 3)!=0;ii++) ;
  if (ii==iflnum) {printf("JitCon ERR: %s not found as a flag (%s)\n",sf,iflags); hxe();}
  if (ifarg(2)) (&ip->type)[ii] = val = (unsigned char) *getarg(2);  
  _lflag=(double)(unsigned char)(&ip->type)[ii];
  if (ifarg(3)) for (i=0;i<cesz;i++) { 
    lop(ce,i); 
    (&qp->type)[ii]=val;
  }
  }
  ENDVERBATIM
}
