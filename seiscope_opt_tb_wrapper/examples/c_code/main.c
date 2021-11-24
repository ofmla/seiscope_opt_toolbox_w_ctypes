#include <stdio.h>
#include <stdlib.h>
#include "optim_bc.h"

int main() {
    int flag;    // communication FLAG
    int n;       // dimension of the problem
    float fcost; // cost function value
    
	printf("Example of C calling Fortran SEISCOPE optimization toolbox subroutines\n");
    
    // declare pointer variables to point to allocated heap space
    float *x;    // current point
    float *grad; // current gradient

	n =2;     // dimension
	flag = 0; // first flag
	
	// call malloc to allocate that appropriate number of bytes for the arrays
	x = (float *)malloc(sizeof(float)*n);  // allocate 2 floats
	grad = (float *)malloc(sizeof(float)*n);  
	
	// Create structure
	optim_type optim;
	optim.niter_max = 10000; // maximum iteration number
	optim.conv = 1e-8;       // tolerance for the stopping criterion
	optim.print_flag = 1;    // print info in output files
	optim.debug = false;     // level of details for output files
	optim.l = 20;            // maximum number of stored pairs used for
                             // the l-BFGS approximation
	
	// initial guess
	x[0] = 1.5;
	x[1] = 1.5;
	
	/* computation of the cost and gradient associated
	 * with the initial guess
	 */
	rosenbrock(x,&fcost,grad);
	
    /* optimization loop: while convergence not reached or 
     * linesearch not failed, iterate
     */
	do
	{   // call fortran
		LBFGS(n,x,&fcost,grad,&optim,&flag,NULL,NULL);
		if (flag == 1) {
			rosenbrock(x,&fcost,grad);}
	}while (flag != 2 && flag != 4);

    // Helpful console writings
    printf("END OF TEST\n");
    printf("FINAL iterate is : %f %f\n", x[0], x[1]);
	if(optim.print_flag){
    printf("See the convergence history in iterate_LB.dat\n");}

    return 0;
}

