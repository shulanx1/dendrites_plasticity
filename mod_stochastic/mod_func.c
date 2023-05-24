#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _K_Pst_2F_reg();
extern void _K_Tst_2F_reg();
extern void _NaTs2_t_2F_reg();
extern void _Nap_Et2_2F_reg();
extern void _na_2F_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," K_Pst_2F.mod");
fprintf(stderr," K_Tst_2F.mod");
fprintf(stderr," NaTs2_t_2F.mod");
fprintf(stderr," Nap_Et2_2F.mod");
fprintf(stderr," na_2F.mod");
fprintf(stderr, "\n");
    }
_K_Pst_2F_reg();
_K_Tst_2F_reg();
_NaTs2_t_2F_reg();
_Nap_Et2_2F_reg();
_na_2F_reg();
}
