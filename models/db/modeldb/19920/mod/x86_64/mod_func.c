#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _ca_a1g_reg(void);
extern void _ca_a1h_reg(void);
extern void _ca_a1i_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," \"ca_a1g.mod\"");
    fprintf(stderr," \"ca_a1h.mod\"");
    fprintf(stderr," \"ca_a1i.mod\"");
    fprintf(stderr, "\n");
  }
  _ca_a1g_reg();
  _ca_a1h_reg();
  _ca_a1i_reg();
}
