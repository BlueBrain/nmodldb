#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _invlfire_reg(void);
extern void _misc_reg(void);
extern void _stats_reg(void);
extern void _vecst_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," \"invlfire.mod\"");
    fprintf(stderr," \"misc.mod\"");
    fprintf(stderr," \"stats.mod\"");
    fprintf(stderr," \"vecst.mod\"");
    fprintf(stderr, "\n");
  }
  _invlfire_reg();
  _misc_reg();
  _stats_reg();
  _vecst_reg();
}
