#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _AMPA_reg(void);
extern void _GABAA_reg(void);
extern void _kdr_reg(void);
extern void _matrix_reg(void);
extern void _naf_reg(void);
extern void _passiv_reg(void);
extern void _precall_reg(void);
extern void _pregen_reg(void);
extern void _pregencv_reg(void);
extern void _pulsecv_reg(void);
extern void _sinstim_reg(void);
extern void _vecst_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," \"AMPA.mod\"");
    fprintf(stderr," \"GABAA.mod\"");
    fprintf(stderr," \"kdr.mod\"");
    fprintf(stderr," \"matrix.mod\"");
    fprintf(stderr," \"naf.mod\"");
    fprintf(stderr," \"passiv.mod\"");
    fprintf(stderr," \"precall.mod\"");
    fprintf(stderr," \"pregen.mod\"");
    fprintf(stderr," \"pregencv.mod\"");
    fprintf(stderr," \"pulsecv.mod\"");
    fprintf(stderr," \"sinstim.mod\"");
    fprintf(stderr," \"vecst.mod\"");
    fprintf(stderr, "\n");
  }
  _AMPA_reg();
  _GABAA_reg();
  _kdr_reg();
  _matrix_reg();
  _naf_reg();
  _passiv_reg();
  _precall_reg();
  _pregen_reg();
  _pregencv_reg();
  _pulsecv_reg();
  _sinstim_reg();
  _vecst_reg();
}
