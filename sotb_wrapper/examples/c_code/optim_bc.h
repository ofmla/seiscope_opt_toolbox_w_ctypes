#include <stdbool.h>
#ifdef __cplusplus
extern "C" {
#endif
#ifndef OPTIM_BC_H
#define OPTIM_BC_H
/* Created a structure here. The name of the structure is
 * optim_type.
 */
typedef struct {
	bool debug;
	float threshold;
	int print_flag;
	bool first_ls;
    int task;
    int nls_max, cpt_ls, nfwd_pb, cpt_iter, niter_max;
    float f0, fk, conv;
    float m1, m2, mult_factor, alpha_L, alpha_R, alpha;
    float q0, q;
    int cpt_lbfgs, l;
    bool conv_CG;
    int cpt_iter_CG, niter_max_CG, nhess;
    int CG_phase;
    int comm;
    float qk_CG, qkm1_CG, hessian_term, eta, norm_grad, norm_grad_m1, norm_residual;
    float dHd, res_scal_respreco, alpha_CG;
} optim_type;
#endif

void rosenbrock(float *, float *, float *);
void rosenbrock_hess(float *, float *, float *);
void set_inputs(optim_type*, float, int, float*, int*, float*, int*, int*, int*, bool*);
void PSTD(int, float*, float*, float*, optim_type*, int*, void*, void*);
void PNLCG(int, float*, float*, float*, optim_type*, int*, void*, void*);
void LBFGS(int, float*, float*, float*, optim_type*, int*, void*, void*);
void PLBFGS(int, float*, float*, float*, float*, float*, optim_type *, int*, void*, void*);
void TRN(int, float*, float*, float*, float*, float*, optim_type*, int*, void*, void*);
void PTRN(int, float*, float*, float*, float*, float*, float*, float*, float*, optim_type*, int*, void*, void*);

#ifdef __cplusplus
}
#endif
