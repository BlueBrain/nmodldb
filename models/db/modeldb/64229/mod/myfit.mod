: $Id: myfit.mod,v 1.2 2006/02/14 13:11:29 hines Exp $
 
:* COMMENT
COMMENT
ENDCOMMENT

NEURON {
  SUFFIX nothing
  GLOBAL FIT_INSTALLED
}

PARAMETER {
  FIT_INSTALLED=0
}

ASSIGNED { RES }

VERBATIM
#include <stdlib.h>
#include <math.h>
/*#include <values.h> /* contains MAXLONG */
#include <limits.h> /* contains MAXLONG */
#include <sys/time.h> 
extern double* hoc_pgetarg();
extern double hoc_call_func(Symbol*, int narg);
extern FILE* hoc_obj_file_arg(int narg);
extern Object** hoc_objgetarg();
extern void vector_resize();
extern int vector_instance_px();
extern void* vector_arg();
extern double* vector_vec();
extern double hoc_epsilon;
extern void set_seed();
extern int ivoc_list_count(Object*);
extern Object* ivoc_list_item(Object*, int);
extern int list_vector_px2();
int list_vector_px();
int list_vector_resize();

typedef struct BVEC {
 int size;
 int bufsize;
 short *x;
 Object* o;
} bvec;
ENDVERBATIM

: spkcmp(modlist,targlist,flag)
: compare spike times from left to right finding the closest for each one
: better algorithm would do all pairwise comparisons and take the smallest, next smallest, etc.
VERBATIM
static double spkcmp (void* vv) {
  int flag, i, j, k, nx, nm[10], nt[10], num, sum, minind;
  Object *ob, *ob1;
  double *vvm[10], *vvt[10], *x, err, diff, min;

  nx = vector_instance_px(vv, &x);
  ob  =   *hoc_objgetarg(1);
  ob1 =   *hoc_objgetarg(2);
  if (ifarg(3)) flag=(int)*getarg(3); else flag=0;

  num = ivoc_list_count(ob);
  i = ivoc_list_count(ob1);
  if (num>10) hoc_execerror("ERR: spkcmp can only handle 10 vectors", 0);
  if (num!=i) hoc_execerror("ERR: spkcmp different sized lists", 0);
  if (nx!=num*2) hoc_execerror("ERR: spkcmp vec should be 2*list length", 0);
  for (i=0;i<num;i++) { 
    nm[i] = list_vector_px(ob,  i,  &vvm[i]);
    nt[i] = list_vector_px(ob1, i,  &vvt[i]);
  }
  for (i=0;i<num;i++) {            /* spike vectors */
    if (nm[i]<=nt[i]) { /* fewer model spikes */
      for (j=0,err=0;j<nm[i];j++) {     /* model spike times */
        for (k=0,min=1e9;k<nt[i];k++) {           /* find min time diff with target spikes */
          diff=fabs(vvm[i][j]-vvt[i][k]);
          if (diff<min) { 
            min=diff;       /* closest spike to this model spike */
            minind=k;
          }
          /* printf("%d %d %d %g %g %g %g %g\n",i,j,k,vvm[i][j],vvt[i][k],diff,min,err); */
        }
        err+=min; 
        vvt[i][minind]=1e9; /* remove that spike */
      }
    } else {  /* go through target spike list */
      for (k=0,err=0;k<nt[i];k++) {         /*  */
        for (j=0,min=1e9;j<nm[i];j++) {     /* model spike times */
          diff=fabs(vvm[i][j]-vvt[i][k]);
          if (diff<min) { 
            min=diff;       /* closest spike to this model spike */
            minind=j;
          }
          /* printf("%d %d %d %g %g %g %g %g\n",i,j,k,vvm[i][j],vvt[i][k],diff,min,err); */
        }
        err+=min; 
        vvm[i][minind]=1e9; /* remove that spike */
      }
    }
    x[i]=err;
  }
  for (i=0;i<num;i++) x[i+num]=(double)abs(nm[i]-nt[i]);
  for (i=0,err=0;i<2*num;i++) err+=x[i];
  return err;
}
ENDVERBATIM

: bursty(modlist,mininvl) -- measure burstiness
VERBATIM
static double bursty (void* vv) {
  int i, j, k, nx, nm[10], nt[10], num, minind;
  Object *ob, *ob1;
  double *vvm[10], *vvt[10], *x, err, diff, mininvl, invl, last, sum, cnt;

  nx = vector_instance_px(vv, &x);
  ob  =   *hoc_objgetarg(1);
  if (ifarg(2)) mininvl=*getarg(2); else mininvl=3.;

  num = ivoc_list_count(ob);
  if (num>10) hoc_execerror("ERR: spkcmp can only handle 10 vectors", 0);
  if (nx!=2*num) hoc_execerror("ERR: spkcmp vec should be 2*list length", 0);
  for (i=0;i<num;i++) { 
    nm[i] = list_vector_px(ob,  i,  &vvm[i]);
  }
  for (i=0;i<num;i++) {            /* spike vectors */
    for (j=1,last=vvm[i][0],sum=0.,cnt=0.;j<nm[i];j++) {
      if ((invl=vvm[i][j]-last)<mininvl) { sum+=invl; cnt+=1; }
      last=vvm[i][j];
    }
    x[2*i]=cnt; x[2*i+1]=sum/cnt; /* #of spikes in bursts, ave invl */
  }
  return x[0];
}
ENDVERBATIM

:* PROCEDURE install_myfit()
PROCEDURE install_myfit () {
  FIT_INSTALLED=1
  VERBATIM
  install_vector_method("spkcmp", spkcmp);
  install_vector_method("bursty", bursty);
  ENDVERBATIM
}

:* idptr(&x) gives location of a pointer
PROCEDURE idptr () {
VERBATIM
  double *fr;
  fr = hoc_pgetarg(1);
  printf("%x\n",fr);
ENDVERBATIM
}

:* veclistcp(src,dest)
FUNCTION veclistcp () {
  VERBATIM {
  int code, i, j, ns, nd, num, n[2];
  Object *obs, *obd;
  double *vvs, *vvd;

  obs =   *hoc_objgetarg(1);
  obd =   *hoc_objgetarg(2);
  num = ivoc_list_count(obs);
  i =   ivoc_list_count(obd);
  if (num!=i) hoc_execerror("veclistcp ERR: different sized lists", 0);

  for (i=0;i<num;i++) { 
    ns = list_vector_px(obs,  i, &vvs);
    nd = list_vector_px(obd, i, &vvd);
    if (ns!=nd) { printf("veclistcp ERR %d %d %d\n",i,ns,nd);
      hoc_execerror("Vectors must all be same size: ", 0); }
    for (j=0;j<ns;j++) vvd[j]=vvs[j];
  }
  return num;
  }
  ENDVERBATIM
}
