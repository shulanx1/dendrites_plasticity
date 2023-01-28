#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _IM_reg();
extern void _exp2syn_NMDA_reg();
extern void _hh2_reg();
extern void _vecevent_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," IM.mod");
fprintf(stderr," exp2syn_NMDA.mod");
fprintf(stderr," hh2.mod");
fprintf(stderr," vecevent.mod");
fprintf(stderr, "\n");
    }
_IM_reg();
_exp2syn_NMDA_reg();
_hh2_reg();
_vecevent_reg();
}
