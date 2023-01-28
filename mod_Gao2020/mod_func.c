#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _CaT_reg();
extern void _Cad_reg();
extern void _IL_reg();
extern void _ca_reg();
extern void _kBK_reg();
extern void _kadist_reg();
extern void _kaprox_reg();
extern void _kv_reg();
extern void _na_reg();
extern void _vmax_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," CaT.mod");
fprintf(stderr," Cad.mod");
fprintf(stderr," IL.mod");
fprintf(stderr," ca.mod");
fprintf(stderr," kBK.mod");
fprintf(stderr," kadist.mod");
fprintf(stderr," kaprox.mod");
fprintf(stderr," kv.mod");
fprintf(stderr," na.mod");
fprintf(stderr," vmax.mod");
fprintf(stderr, "\n");
    }
_CaT_reg();
_Cad_reg();
_IL_reg();
_ca_reg();
_kBK_reg();
_kadist_reg();
_kaprox_reg();
_kv_reg();
_na_reg();
_vmax_reg();
}
