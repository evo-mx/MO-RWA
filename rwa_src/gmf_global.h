/*
 * gmf_global.h
 *
 *  Created on: Jul 22, 2020
 *      Author: Saul Zapotecas
 */
#ifndef GMF_GLOBAL_H_
#define GMF_GLOBAL_H_

#define gmf_true 1
#define gmf_false 0

extern void (*gmf_test_problem)(double *F, double *G, double *xr, int *xi,
		int *xb);


struct
{
	char benchmarck[50];
	char name[50];

	char PF_file[255];
	char PS_file[255];

	size_t nreal;
	size_t nint;
	size_t nbin;

	size_t ncons;
	size_t nobjs;
	size_t dynamic_flag;

	size_t evaluations;
	size_t std;

	double *xmin_real;
	double *xmax_real;
	int *xmin_int;
	int *xmax_int;
} gmf_mop;

#endif /* GMF_GLOBAL_H_ */
