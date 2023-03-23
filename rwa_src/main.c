/*
 * main.c
 *
 *  Created on: Jun 7, 2022
 *      Author: Saul Zapotecas
 */
#include <stdio.h>
#include <stdlib.h>
#include "gmf_global.h"
#include "gmf_dv_rwa.h"

void (*gmf_test_problem)(double *F, double *G, double *xr, int *xi, int *xb);

/**
 * Display decision variables and objectives
 */
void display_solution(double *x, double *f)
{
	int i;
	printf("x: ");
	for (i = 0; i < gmf_mop.nreal; ++i)
	{
		printf("%lf ", x[i]);
	}
	printf("\nf: ");
	for (i = 0; i < gmf_mop.nobjs; ++i)
	{
		printf("%lf ", f[i]);
	}
	printf("\n");
	return;
}

/**
 * Main Function
 */
int main()
{
	double *x, *f;

	srand(10000);

	/* Usage example. Function Subasi2016 */
	gmf_rwa_setup("Subasi2016");
	printf("Function: Subasi2006 (%zu vars. y %zu objs.)", gmf_mop.nreal,
			gmf_mop.nobjs);
	x = malloc(sizeof(double) * gmf_mop.nreal);
	f = malloc(sizeof(double) * gmf_mop.nobjs);
	gmf_rwa_rnd_solution(x);
	gmf_test_problem(f, NULL, x, NULL, NULL);
	display_solution(x, f);
	gmf_rwa_setdown();
	free(x);
	free(f);
	printf("\n");

	/* Usage example. Function Ahmad2017 */
	gmf_rwa_setup("Ahmad2017");
	printf("Function: Ahmad2017 (%zu vars. y %zu objs.)", gmf_mop.nreal,
			gmf_mop.nobjs);
	x = malloc(sizeof(double) * gmf_mop.nreal);
	f = malloc(sizeof(double) * gmf_mop.nobjs);
	gmf_rwa_rnd_solution(x);
	gmf_test_problem(f, NULL, x, NULL, NULL);
	display_solution(x, f);
	gmf_rwa_setdown();
	free(x);
	free(f);
	printf("\n");

	return 0;
}

