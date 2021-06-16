: $Id: bnet.mod,v 1.49 2012/08/09 16:26:55 samn Exp $ 

:* main COMMENT
COMMENT

boolean network 

ENDCOMMENT

:* main VERBATIM block
VERBATIM

#include "misc.h"

// BP gets a pointer to the boonet* that belongs to each BNET object
#define BP (*((boonet**) &(_p_sop)))

// node (can represent a molecule)
typedef struct BNODE {
  char* name;    // name of the node
  int id;        // id
  int state;  // state of the node
  int count;  // temporary counter
  int knockout; // is this node knocked-out? 
  int start; // starting state
  int sthresh; // special threshold for switching the node's state -- 0 means it's a normal node
  int scount; // special count (# of times its rule is activated)
  double x,y,z;  // location of the node
} bnode;

// represents a rule (contains edges or links between nodes)
typedef struct BRULE {
  struct BNODE** psrc; // source (usually 1 node but can be ANDed)
  struct BNODE*  ptarg;   // pointer to target
  int weight;             // weight of the rule ( < 0 == inhib, > 0 == activate )
  int* psrcstate;       // state of source nodes required for this rule
  int nsrc;           // # of sources in the rule
} brule;

// C-level structure/representation of the boolean network (stores nodes)
typedef struct BOONET {
  struct BNODE* pnodes;
  int numnodes;
  int id;
  struct BRULE* prules; // list of outgoing rules (pointers to targets and weights)
  int nrules;  // number of outgoing edges
  int rulebufsz; // prules buffer size (in units of brules)
} boonet;

// set the state of each node to its starting value
void startboonet (boonet* pnet) {
  int i;
  if(verbose) printf("startboonet\n");
  for(i=0;i<pnet->numnodes;i++) {
    pnet->pnodes[i].count = pnet->pnodes[i].scount = 0;
    pnet->pnodes[i].state = pnet->pnodes[i].start; // set to starting
  }
}

// create and return the C-level structure/representation of the boolean network
boonet* makeboonet (int numnodes) {
  int i;
  boonet* pnet;
  pnet = (boonet*) calloc(1, sizeof(boonet)); // allocate memory for the boonet
  pnet->pnodes = (bnode*) calloc(numnodes, sizeof(bnode)); // allocate memory for the nodes
  pnet->numnodes = numnodes; // store # of nodes
  startboonet(pnet);         // set states to starting values
  for(i=0;i<numnodes;i++) pnet->pnodes[i].id = i; // assign the IDs
  return pnet;
}

// free a single bnode
void freebnode(bnode* p) {
  if(p->name) {
    free(p->name);
    p->name = 0x0;
  }
}

// free a single brule's memory
void freebrule(brule* p) {
  if(p->psrc) free(p->psrc);
  if(p->psrcstate) free(p->psrcstate);
}

// free the entire boolean network
void freeboonet (boonet* pnet) {
  int i;
  if(pnet->pnodes) {
    for(i=0;i<pnet->numnodes;i++) freebnode(&pnet->pnodes[i]);
    free(pnet->pnodes);
  }  
  pnet->pnodes = 0x0;
  if(pnet->prules) {
    for(i=0;i<pnet->nrules;i++) freebrule(&pnet->prules[i]);
    free(pnet->prules);
    pnet->prules = 0x0;
  }
  free(pnet);
}

// add a rule from psrc to targid with specified weight and source state
// return pointer to the rule
brule* addrule (boonet* pnet, double* psrc, int targid, double weight, double* psrcstate, int nsrc) {
  int idx,i;
  if(pnet->prules == 0x0) { // take care of memory allocation
    pnet->rulebufsz = 16;
    pnet->prules = calloc(pnet->rulebufsz, sizeof(brule) );
    pnet->nrules = 0;
  } else if(pnet->nrules >= pnet->rulebufsz) {
    pnet->rulebufsz *= 2;
    pnet->prules = realloc(pnet->prules, pnet->rulebufsz * sizeof(brule) );
  }
  idx = pnet->nrules;
  pnet->prules[idx].nsrc = nsrc;           // set # of sources
  pnet->prules[idx].ptarg = &pnet->pnodes[targid]; // set target
  pnet->prules[idx].weight = weight;               // set weight
  // set the sources
  pnet->prules[idx].psrc = calloc(nsrc, sizeof(bnode*));
  for(i=0;i<nsrc;i++) pnet->prules[idx].psrc[i] = &pnet->pnodes[(int)psrc[i]];
  // set source state for rule to be 'on'
  pnet->prules[idx].psrcstate = calloc(nsrc, sizeof(int));
  for(i=0;i<nsrc;i++) pnet->prules[idx].psrcstate[i] = (int) psrcstate[i]; 
  return &pnet->prules[ pnet->nrules++ ]; // return the new rule and inc # rules
}

// advance the full boolean network by a single iteration
void advanceboonet (boonet* pnet) {
  int i , j, nn = pnet->numnodes, nsrcact;
  bnode *pnodes = pnet->pnodes;
  brule *prule;
  for(i=0;i<nn;i++) pnodes[i].count = 0;      // initialize temporary counts to 0
  for(i=0;i<pnet->nrules;i++) {
    prule = &pnet->prules[i];
    nsrcact = 0; // number of sources that are activated
    for(j=0;j<prule->nsrc;j++) { // go thru the sources of the rule
      if(prule->psrc[j]->knockout) continue; // skip knockouts
      if(prule->psrc[j]->state == prule->psrcstate[j]){//check if source node in correct state
        nsrcact += 1; // count number of source nodes in correct state (for activation of the rule)
      }
    }
    if(prule->nsrc == nsrcact) prule->ptarg->count += prule->weight;//integrate if ALL sources in proper state
  }
  for(i=0;i<nn;i++) { // update the state of each node based on its inputs
    if(pnodes[i].knockout) {pnodes[i].scount=pnodes[i].count=pnodes[i].state=0; continue;}
    if(pnodes[i].count > 0) {       // turn on (more activating vs inhibiting input)     
      if(pnodes[i].sthresh > 0) { // special node?
        pnodes[i].scount += 1;    // increase the special counter
      } else {//normal node so gets turned on when activation count > inhibition count
        pnodes[i].state = 1;
      }      
    } else if(pnodes[i].count<0||(pnodes[i].sthresh>0 && pnodes[i].count==0)){//turn off (inhibitory > activating input) or special node
      pnodes[i].state = 0;
      pnodes[i].scount = 0; // special node count reset to 0
    }     
    if(pnodes[i].sthresh > 0){//set the state of any special nodes here
      if(pnodes[i].scount >= pnodes[i].sthresh) { // passed the threshold? turn on
        pnodes[i].state = 1;
      } else pnodes[i].state = 0;
    }
  }
}

ENDVERBATIM

:* NEURON, PARAMETER, ASSIGNED blocks
NEURON {
  ARTIFICIAL_CELL BNET
  : POINT_PROCESS BNET
  RANGE tstep
  RANGE xloc, yloc, zloc
  GLOBAL verbose, installed
  POINTER sop :::: Structure pointer for the boolean network
}

PARAMETER {
  tstep = 0
  verbose = 0
  sop = 0
}

ASSIGNED {
  installed
}

:* CONSTRUCTOR, DESTRUCTOR, INITIAL
:** CONSTRUCT: create a structure to save the identity of this unit and char integer flags
CONSTRUCTOR {
  VERBATIM 
  boonet* pnet;
  int sz; 
  if((sz = (int)*getarg(1)) < 1) {
    printf("BNET err0: must have a network with positive # of nodes!\n");
    hxe();
  }
  _p_sop = (void*) makeboonet(sz);
  pnet = BP;
  pnet->id = ifarg(2) ? (int) *getarg(2) : 0;
  ENDVERBATIM
}

: setrule(vsource, targid, weight, vsourcestate) - sets the rule from a set of 3 vectors
: vsource has node IDs of sources, targid is ID of target, weight is weight
: vsourcestate has states the source must be in for the given rule to be turned on.
FUNCTION setrule () {
  VERBATIM
  double *psrc, weight, *psrcstate;
  int targid, i, nsrc;
  if (! ifarg(4) ) {
    printf("BNET.setrule(vsource, targid, weight, vsourcestates) - sets a rule.\n");
    printf("vsource has node IDs of sources, targid is ID of target, weight is (-1,1) for inhib/activating\n");
    printf("vsourcestate has states the source must be in for the given rule to be turned on.\n");
    return 0.0;
  }
  if ( (nsrc = vector_arg_px(1,&psrc)) < 1) {
    printf("BNET.setrule WARN0: empty source Vector!\n");
    return 0.0;
  }
  targid = (int) *getarg(2);
  if(targid < 0 || targid >= BP->numnodes) {
    printf("BNET.setrule ERR0: invalid target id : %d\n",targid);
    return 0.0;
  }
  weight = (int) *getarg(3);
  if( nsrc != vector_arg_px(4,&psrcstate) ) {
    printf("BNET.setrule ERR1: vsource, vsourcestate must have same size!\n");
    return 0.0;
  }
  for(i=0;i<nsrc;i++) {
    if( psrc[i] < 0 || psrc[i] >= BP->numnodes) {
      printf("BNET.setrule ERR2: invalid source node id %d. netsize=%d\n",(int)psrc[i],BP->numnodes);
      return 0.0;
    }
    if(verbose>1) printf("adding rule from %d -> %d : w = %g\n",(int)psrc[i],targid,weight);
  }
  addrule(BP, psrc, targid, weight, psrcstate, nsrc);
  return 1.0;
  ENDVERBATIM
}

: get rid of all the edges (but keep the nodes)
PROCEDURE clearrules () {
  VERBATIM
  BP->nrules = 0;
  ENDVERBATIM
}

: print the network
PROCEDURE pr () {
  VERBATIM
  int i, j, k;
  bnode* pnodes = BP->pnodes;
  brule* prule;
  char srcstr[4096], stmp[4096];
  printf("net: numnodes=%d, numrules=%d\n",BP->numnodes,BP->nrules);
  for(i=0;i<BP->numnodes;i++) {
    if(pnodes[i].name) {
      printf("%s: state=%d, count=%d\n",pnodes[i].name,pnodes[i].state,pnodes[i].count);
    } else {
      printf("%d: state=%d, count=%d\n",i,pnodes[i].state,pnodes[i].count);
    }    
  }
  for(i=0;i<BP->nrules;i++) {
    prule = &BP->prules[i];
    srcstr[0]=0;
    for(j=0;j<prule->nsrc;j++) {
      if(prule->psrc[j]->name) {
        sprintf(stmp,"%s%s%s ", j>0?"AND ":"", prule->psrcstate[j]?"":"!", prule->psrc[j]->name);
      } else {
        sprintf(stmp,"%s%s%s ", j>0?"AND ":"", prule->psrcstate[j]?"":"!", prule->psrc[j]->id);
      }
      strcat(srcstr,stmp);
    }
    if(prule->ptarg->name) {
      printf("%s-> %s , w = %d\n",srcstr,prule->ptarg->name,prule->weight);
    } else {
      printf("%s-> %d , w = %d]\n",srcstr,prule->ptarg->id,prule->weight);
    }
  }
  ENDVERBATIM
}

: BNET.graphviz([dotname,imageename,format,L->R direction,width,height,fontsize]) - print the network as a graphviz string,
: and optionally save to a dot file (dotname) and to png/pdf. format should be png or pdf or other output formats
: supported by graphviz.
FUNCTION graphviz () {
  VERBATIM
  int i, j, k, LR, fsz, w, h;
  bnode* pnodes = BP->pnodes;
  brule* prule;
  char *ncolor, *fcolor, *arrowtype, *lstyle, *shape;//node color, font color, arrow type, line style, node shape
  char buf[4096], *dotname, *fname, *ext, fontsize[128];
  double penw; // penwidth
  FILE* fp = 0x0;
  dotname = ifarg(1) ? gargstr(1) : 0x0;
  fname = ifarg(2) ? gargstr(2) : 0x0;
  ext =   ifarg(3) ? gargstr(3) : "gif";
  LR = ifarg(4) ? (int) *getarg(4) : 1;
  w = ifarg(5) ? (int) *getarg(5) : -1;
  h = ifarg(6) ? (int) *getarg(6) : -1;
  fsz = ifarg(7) ? (int) *getarg(7) : -1;
  if(fsz==-1) sprintf(fontsize,"%s"," "); else sprintf(fontsize,"fontsize=%d,",fsz);
  if(fname) if( !(fp = fopen(dotname,"w"))) {
    printf("BNET.graphviz ERR0: could not open %s\n",fname);
    return 0.0;
  }
  sprintf(buf, "%s", "digraph G {\n"); if(fp) fprintf(fp,"%s",buf); else fprintf(stdout,"%s",buf); 
  if(LR){sprintf(buf, "%s", "\trankdir=LR;\n"); if(fp) fprintf(fp,"%s",buf); else fprintf(stdout,"%s",buf); }
  if(w>0 && h>0) {sprintf(buf, "size=\"%d,%d\"\n",w,h); if(fp) fprintf(fp,"%s",buf); else fprintf(stdout,"%s",buf);}
  for(i=0;i<BP->numnodes;i++) {
    ncolor = BP->pnodes[i].knockout ? "white" : BP->pnodes[i].state > 0 ? "black" : "gray";
    fcolor = BP->pnodes[i].knockout ? "black" : "white";
    shape = pnodes[i].sthresh > 0 ? "invtriangle" : "doublecircle";
    if(BP->pnodes[i].name) {
      sprintf(buf,"\t%s [fontcolor=%s,%sstyle=filled,shape=%s,fillcolor=%s,color=%s]\n",
              BP->pnodes[i].name,fcolor,fontsize,shape,ncolor,ncolor);
    } else {
      sprintf(buf,"\t%d [fontcolor=%s,%sstyle=filled,shape=%s,fillcolor=%s,color=%s]\n",
              i,fcolor,fontsize,shape,ncolor,ncolor);
    }
    if(fp) fprintf(fp,"%s",buf); else fprintf(stdout,"%s",buf); 
  }
  for(i=0;i<BP->nrules;i++) {
    prule = &BP->prules[i];
    for(j=0;j<prule->nsrc;j++) {
      penw = prule->psrcstate[j] == prule->psrc[j]->state ? 6.0 : 1.0;
      arrowtype = prule->weight < 0 ? "tee" : "open";
      lstyle = prule->psrcstate[j] == 0 ? ",style=dashed" : " ";
      if(prule->psrc[j]->name) {
        if(prule->ptarg->name) {
          sprintf(buf,"\t%s -> %s [arrowhead=%s,penwidth=%g,color=%s%s]\n",
                  prule->psrc[j]->name,prule->ptarg->name,arrowtype,penw,prule->weight>0?"red":"blue",lstyle);
        } else {
          sprintf(buf,"\t%s -> %d [arrowhead=%s,penwidth=%g,color=%s%s]\n",
                  prule->psrc[j]->name,prule->ptarg->id,arrowtype,penw,prule->weight>0?"red":"blue",lstyle);
        }
      } else if(prule->ptarg->name) {
          sprintf(buf,"\t%d -> %s [arrowhead=%s,penwidth=%g,color=%s%s]\n",
                  prule->psrc[j]->id,prule->ptarg->name,arrowtype,penw,prule->weight>0?"red":"blue",lstyle);
      } else {
          sprintf(buf,"\t%d -> %d [arrowhead=%s,penwidth=%g,color=%s%s]\n",
                  prule->psrc[j]->id,prule->ptarg->id,arrowtype,penw,prule->weight>0?"red":"blue",lstyle);
      }
      if(fp) fprintf(fp,"%s",buf); else fprintf(stdout,"%s",buf); 
    }
  }
  sprintf(buf,"%s","}\n"); if(fp) fprintf(fp,"%s",buf); else fprintf(stdout,"%s",buf); 
  if(fp) fclose(fp);
  if(fname) {
    sprintf(buf,"dot %s -T%s > %s",dotname,ext,fname);
    if(0!=system(buf)) {printf("BNET.graphviz ERR1 : couldn't run %s\n",buf); return 0.0;}
  }
  return 1.0;
  ENDVERBATIM
}

DESTRUCTOR {
  VERBATIM
  freeboonet(BP);
  ENDVERBATIM
}

: all nodes in network set to starting states
PROCEDURE start () {
  tstep = 0
  VERBATIM
  startboonet(BP); 
  ENDVERBATIM
}

:** INITIAL
INITIAL {
  start()
}

: strvalid - return the index to be used for get/setnodevals
FUNCTION strvalid () {
  VERBATIM
  char *pname;
  static char *pnames[6] = {"state", "count", "knockout", "start", "sthresh", "scount"};
  int i;
  pname = gargstr(1);
  for(i=0;i<6;i++) {
    if(!strcmp(pname,pnames[i])) return i;
  }
  return -1;
  ENDVERBATIM
}

:** getscount(vec) - retrieves BNET scount into a vector. each element indexed into the vector
: corresponds to the node with the given id.
PROCEDURE getscount () {
  VERBATIM
  double *ps; int i; void *vs;
  if(!ifarg(1)) {
    printf("BNET.getscount(vec) - returns scount of each node in vec\n");
    return 0; 
  }
  vs = vector_arg(1);
  ps = vector_newsize(vs,BP->numnodes);
  for(i=0;i<BP->numnodes;i++) ps[i] = (double) BP->pnodes[i].scount;
  ENDVERBATIM
}

:** setscount(vec) - sets vector vec into BNET node scount. each element indexed into the vector
: corresponds to the node with the given id.
FUNCTION setscount () {
  VERBATIM
  double *ps; int i, sz;
  if(!ifarg(1)) {
    printf("BNET.setscount(vec) - sets scount of each node in vec\n");
    return 0.0;
  }
  if( (sz = vector_arg_px(1,&ps)) != BP->numnodes ) {
    printf("BNET.setscount ERR0: vec.size(%d) != BNET.numnodes(%d)\n",sz,BP->numnodes);
    return 0.0;
  }
  for(i=0;i<BP->numnodes;i++) BP->pnodes[i].scount = (int) ps[i];
  return 1.0;
  ENDVERBATIM
}

:** getsthresh(vec) - retrieves BNET sthresh into a vector. each element indexed into the vector
: corresponds to the node with the given id.
PROCEDURE getsthresh () {
  VERBATIM
  double *ps; int i; void *vs;
  if(!ifarg(1)) {
    printf("BNET.getsthresh(vec) - returns sthresh of each node in vec\n");
    return 0;
  }
  vs = vector_arg(1);
  ps = vector_newsize(vs,BP->numnodes);
  for(i=0;i<BP->numnodes;i++) ps[i] = (double) BP->pnodes[i].sthresh;
  ENDVERBATIM
}

:** setsthresh(vec) - sets vector vec into BNET node sthresh. each element indexed into the vector
: corresponds to the node with the given id.
FUNCTION setsthresh () {
  VERBATIM
  double *ps; int i, sz;
  if(!ifarg(1)) {
    printf("BNET.setsthresh(vec) - sets sthresh of each node in vec\n");
    return 0.0;
  }
  if( (sz = vector_arg_px(1,&ps)) != BP->numnodes ) {
    printf("BNET.setsthresh ERR0: vec.size(%d) != BNET.numnodes(%d)\n",sz,BP->numnodes);
    return 0.0;
  }
  for(i=0;i<BP->numnodes;i++) BP->pnodes[i].sthresh = (int) ps[i];
  return 1.0;
  ENDVERBATIM
}

:** getcount(vec) - retrieves BNET count into a vector. each element indexed into the vector
: corresponds to the node with the given id.
PROCEDURE getcount () {
  VERBATIM
  double *ps; int i; void *vs;
  if(!ifarg(1)) {
    printf("BNET.getcount(vec) - returns count of each node in vec\n");
    return 0;
  }
  vs = vector_arg(1);
  ps = vector_newsize(vs,BP->numnodes);
  for(i=0;i<BP->numnodes;i++) ps[i] = (double) BP->pnodes[i].count;
  ENDVERBATIM
}

:** setcount(vec) - sets vector vec into BNET node counts. each element indexed into the vector
: corresponds to the node with the given id.
FUNCTION setcount () {
  VERBATIM
  double *ps; int i, sz;
  if(!ifarg(1)) {
    printf("BNET.setcount(vec) - sets count of each node in vec\n");
    return 0.0;
  }
  if( (sz = vector_arg_px(1,&ps)) != BP->numnodes ) {
    printf("BNET.setcount ERR0: vec.size(%d) != BNET.numnodes(%d)\n",sz,BP->numnodes);
    return 0.0;
  }
  for(i=0;i<BP->numnodes;i++) BP->pnodes[i].count = (int) ps[i];
  return 1.0;
  ENDVERBATIM
}

:** getstate(vec) - retrieves BNET state into a vector. each element indexed into the vector
: corresponds to the node with the given id.
PROCEDURE getstate () {
  VERBATIM
  double *ps; int i; void *vs;
  if(!ifarg(1)) {
    printf("BNET.getstate(vec) - returns state of each node in vec\n");
    return 0;
  }
  vs = vector_arg(1);
  ps = vector_newsize(vs,BP->numnodes);
  for(i=0;i<BP->numnodes;i++) ps[i] = (double) BP->pnodes[i].state;
  ENDVERBATIM
}

:** setstate(vec) - sets vector vec into BNET node states. each element indexed into the vector
: corresponds to the node with the given id.
FUNCTION setstate () {
  VERBATIM
  double *ps; int i, sz;
  if(!ifarg(1)) {
    printf("BNET.setstate(vec) - sets state of each node in vec\n");
    return 0.0;
  }
  if( (sz = vector_arg_px(1,&ps)) != BP->numnodes ) {
    printf("BNET.setstate ERR0: vec.size(%d) != BNET.numnodes(%d)\n",sz,BP->numnodes);
    return 0.0;
  }
  for(i=0;i<BP->numnodes;i++) BP->pnodes[i].state = (int) ps[i];
  return 1.0;
  ENDVERBATIM
}

:** getstart(vec) - retrieves BNET start states into a vector. each element indexed into the vector
: corresponds to the node with the given id.
PROCEDURE getstart () {
  VERBATIM
  double *ps; int i; void *vs;
  if(!ifarg(1)) {
    printf("BNET.getstart(vec) - returns start state of each node in vec\n");
    return 0;
  }
  vs = vector_arg(1);
  ps = vector_newsize(vs,BP->numnodes);
  for(i=0;i<BP->numnodes;i++) ps[i] = (double) BP->pnodes[i].start;
  ENDVERBATIM
}

:** setstart(vec) - sets vector vec into BNET node start states. each element indexed into the vector
: corresponds to the node with the given id.
FUNCTION setstart () {
  VERBATIM
  double *ps; int i, sz;
  if(!ifarg(1)) {
    printf("BNET.setstart(vec) - sets start state of each node in vec\n");
    return 0.0;
  }
  if( (sz = vector_arg_px(1,&ps)) != BP->numnodes ) {
    printf("BNET.setstart ERR0: vec.size(%d) != BNET.numnodes(%d)\n",sz,BP->numnodes);
    return 0.0;
  }
  for(i=0;i<BP->numnodes;i++) BP->pnodes[i].start = (int) ps[i];
  return 1.0;
  ENDVERBATIM
}

:** setknockout(vec) - sets vector vec into BNET node knockout flag variables. each element indexed into the vector
: corresponds to the node with the given id.
FUNCTION setknockout () {
  VERBATIM
  double *pk; int i, sz;
  if(!ifarg(1)) {
    printf("BNET.knockout(vec) - sets knockout flag of each node in vec\n");
    return 0.0;
  }
  if( (sz = vector_arg_px(1,&pk)) != BP->numnodes ) {
    printf("BNET.knockout ERR0: vec.size(%d) != BNET.numnodes(%d)\n",sz,BP->numnodes);
    return 0.0;
  }
  for(i=0;i<BP->numnodes;i++) BP->pnodes[i].knockout = (int) pk[i];
  return 1.0;
  ENDVERBATIM
}

:** getknockout(vec) - retrieves BNET knockout flags into a vector. each element indexed into the vector
: corresponds to the node with the given id.
PROCEDURE getknockout () {
  VERBATIM
  double *pk; int i; void *vk;
  if(!ifarg(1)) {
    printf("BNET.getknockout(vec) - returns knockout flag of each node in vec\n");
    return 0;
  }
  vk = vector_arg(1);
  pk = vector_newsize(vk,BP->numnodes);
  for(i=0;i<BP->numnodes;i++) pk[i] = (double) BP->pnodes[i].knockout;
  ENDVERBATIM
}

: BNET.setnname(node id, string name) - set node name
FUNCTION setnname () {
  VERBATIM
  int id, sz; char *name;
  id = (int) *getarg(1);
  if(id < 0 || id >= BP->numnodes) {
    printf("BNET.setnname ERR0: invalid node index %d\n",id);
    return 0.0;    
  }
  name = gargstr(2);
  if(!(sz=strlen(name))) {
    printf("BNET.setnname ERR1: empty string\n");
    return 0.0;
  }
  if(BP->pnodes[id].name) free(BP->pnodes[id].name);
  if(!(BP->pnodes[id].name = (char*) malloc(sizeof(char) * (sz + 1)))) {
    printf("BNET.setnname ERR2: couldn't alloc mem for %s\n",name);
    return 0.0;
  }
  strcpy(BP->pnodes[id].name,name);
  return 1.0;
  ENDVERBATIM
}

: getnname(node id, string) - get node name
FUNCTION getnname () {
  VERBATIM
  int i, id, sz; char **pname, string[BUFSIZ];
  char** hoc_pgargstr();
  id = (int) *getarg(1);
  if(id < 0 || id >= BP->numnodes) {
    printf("BNET.getnname ERR0: invalid node index %d\n",id);
    return 0.0;    
  }
  if(!BP->pnodes[id].name || !(sz=strlen(BP->pnodes[id].name))) {
    printf("BNET.getnname ERR1: node %d has no name\n",id);
    return 0.0;
  }
  for(i=0;i<sz && i<BUFSIZ;i++) string[i] = BP->pnodes[id].name[i];
  if(i < BUFSIZ) string[i]=0;
  printf("Aname is %s, %s\n",BP->pnodes[id].name,string);
  pname = hoc_pgargstr(2);
  printf("Bname is %s, %s\n",BP->pnodes[id].name,string);
  hoc_assign_str(pname,string);
  return 1.0;
  ENDVERBATIM
}

: advance the BNET by a single iteration - called from a simulation loop
FUNCTION advancebn () {
  VERBATIM
  advanceboonet(BP);  
  tstep = tstep + 1;
  return tstep;
  ENDVERBATIM
}

FUNCTION advancebnfor () {
  VERBATIM
  int i, n;
  n = (int) *getarg(1);
  for(i=0;i<n;i++) advancebn();
  return tstep;
  ENDVERBATIM
}

FUNCTION numnodes () {
  VERBATIM
  return BP->numnodes;
  ENDVERBATIM
}

FUNCTION numrules () {
  VERBATIM
  return BP->nrules;
  ENDVERBATIM
}

: identifier for BNET
FUNCTION id () {
  VERBATIM
  return (double) BP->id;
  ENDVERBATIM
}

