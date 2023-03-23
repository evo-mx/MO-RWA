/*
 * gmf_dv_rwa.c
 *
 *  Created on: Jun 7, 2022
 *      Author: saul
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "gmf_global.h"

/**
 * According to the description in ref (Subasi et al., 2016)
 * "Multi-objective optimization of a honeycomb heat sink using Response Surface Method"
 * 2 objectives
 * 5 variables
 */
void Subasi2016_setup()
{
	gmf_mop.nobjs = 2;
	gmf_mop.nreal = 5;
	gmf_mop.xmin_real = malloc(sizeof(double) * gmf_mop.nreal);
	gmf_mop.xmax_real = malloc(sizeof(double) * gmf_mop.nreal);
	gmf_mop.xmin_real[0] = 20.0;
	gmf_mop.xmin_real[1] = 6.0;
	gmf_mop.xmin_real[2] = 20.0;
	gmf_mop.xmin_real[3] = 0.0;
	gmf_mop.xmin_real[4] = 8000.0;

	gmf_mop.xmax_real[0] = 60.0;
	gmf_mop.xmax_real[1] = 15.0;
	gmf_mop.xmax_real[2] = 40.0;
	gmf_mop.xmax_real[3] = 30.0;
	gmf_mop.xmax_real[4] = 25000.0;
	return;
}
void Subasi2016(double *F, double *G, double *xr, int *xi, int *xb)
{
	double H = xr[0];
	double t = xr[1];
	double Sy = xr[2];
	double theta = xr[3];
	double Re = xr[4];

	double f;
	double Nu;

	Nu = 89.027 + 0.300 * H - 0.096 * t - 1.124 * Sy - 0.968 * theta
			+ 4.148 * 10e-3 * Re + 0.0464 * H * t - 0.0244 * H * Sy
			+ 0.0159 * H * theta + 4.151 * 10e-5 * H * Re + 0.1111 * t * Sy
			- 4.121 * 10e-5 * Sy * Re + 4.192 * 10e-5 * theta * Re;

	f = 0.4753 - 0.0181 * H + 0.0420 * t + 5.481 * 10e-3 * Sy - 0.0191 * theta
			- 3.416 * 10e-6 * Re - 8.851 * 10e-4 * H * Sy
			+ 8.702 * 10e-4 * H * theta + 1.536 * 10e-3 * t * theta
			- 2.761 * 10e-6 * t * Re - 4.400 * 10e-4 * Sy * theta
			+ 9.714 * 10e-7 * Sy * Re + 6.777 * 10e-4 * H * H;
	F[0] = -Nu;
	F[1] = f;
}

/**
 * According to the description in ref (Goel et al., 2007)
 * "Response surface approximation of Pareto optimal front in multi-objective optimization"
 * 3 objectives
 * 4 variables
 */
void Goel2007_setup()
{
	int i;
	gmf_mop.nobjs = 3;
	gmf_mop.nreal = 4;
	gmf_mop.xmin_real = malloc(sizeof(double) * gmf_mop.nreal);
	gmf_mop.xmax_real = malloc(sizeof(double) * gmf_mop.nreal);
	for (i = 0; i < gmf_mop.nreal; ++i)
	{
		gmf_mop.xmin_real[i] = 0.0;
		gmf_mop.xmax_real[i] = 1.0;
	}
	return;
}

void Goel2007(double *F, double *G, double *xr, int *xi, int *xb)
{
	double a = xr[0];
	double DHA = xr[1];
	double DOA = xr[2];
	double OPTT = xr[3];
	double Xcc;
	double TFmax;
	double TTmax;
	/* Xcc (minimization) */
	Xcc = 0.153 - 0.322 * a + 0.396 * DHA + 0.424 * DOA + 0.0226 * OPTT
			+ 0.175 * a * a + 0.0185 * DHA * a - 0.0701 * DHA * DHA
			- 0.251 * DOA * a + 0.179 * DOA * DHA + 0.0150 * DOA * DOA
			+ 0.0134 * OPTT * a + 0.0296 * OPTT * DHA + 0.0752 * OPTT * DOA
			+ 0.0192 * OPTT * OPTT;
	/* TFmax (minimization) */
	TFmax = 0.692 + 0.477 * a - 0.687 * DHA - 0.080 * DOA - 0.0650 * OPTT
			- 0.167 * a * a - 0.0129 * DHA * a + 0.0796 * DHA * DHA
			- 0.0634 * DOA * a - 0.0257 * DOA * DHA + 0.0877 * DOA * DOA
			- 0.0521 * OPTT * a + 0.00156 * OPTT * DHA + 0.00198 * OPTT * DOA
			+ 0.0184 * OPTT * OPTT;
	/* TTmax (minimization) */
	TTmax = 0.370 - 0.205 * a + 0.0307 * DHA + 0.108 * DOA + 1.019 * OPTT
			- 0.135 * a * a + 0.0141 * DHA * a + 0.0998 * DHA * DHA
			+ 0.208 * DOA * a - 0.0301 * DOA * DHA - 0.226 * DOA * DOA
			+ 0.353 * OPTT * a - 0.0497 * OPTT * DOA - 0.423 * OPTT * OPTT
			+ 0.202 * DHA * a * a - 0.281 * DOA * a * a - 0.342 * DHA * DHA * a
			- 0.245 * DHA * DHA * DOA + 0.281 * DOA * DOA * DHA
			- 0.184 * OPTT * OPTT * a - 0.281 * DHA * a * DOA;

	F[0] = Xcc;
	F[1] = TFmax;
	F[2] = TTmax;
	return;
}

/**
 * According to the description in ref (Liao et al., 2008)
 * "Multiobjective optimization for crash safety design of vehicles using stepwise regression model"
 * 3 objectives
 * 5 variables
 */
void Liao2008_setup()
{
	int i;
	gmf_mop.nobjs = 3;
	gmf_mop.nreal = 5;
	gmf_mop.xmin_real = malloc(sizeof(double) * gmf_mop.nreal);
	gmf_mop.xmax_real = malloc(sizeof(double) * gmf_mop.nreal);
	for (i = 0; i < gmf_mop.nreal; ++i)
	{
		gmf_mop.xmin_real[i] = 1.0;
		gmf_mop.xmax_real[i] = 3.0;
	}
	return;
}
void Liao2008(double *F, double *G, double *xr, int *xi, int *xb)
{
	double Mass;
	double Ain;
	double Intrusion;
	double t1 = xr[0];
	double t2 = xr[1];
	double t3 = xr[2];
	double t4 = xr[3];
	double t5 = xr[4];
	Mass = 1640.2823 + 2.3573285 * t1 + 2.3220035 * t2 + 4.5688768 * t3
			+ 7.7213633 * t4 + 4.4559504 * t5;
	Ain = 6.5856 + 1.15 * t1 - 1.0427 * t2 + 0.9738 * t3 + 0.8364 * t4
			- 0.3695 * t1 * t4 + 0.0861 * t1 * t5 + 0.3628 * t2 * t4
			- 0.1106 * t1 * t1 - 0.3437 * t3 * t3 + 0.1764 * t4 * t4;
	Intrusion = -0.0551 + 0.0181 * t1 + 0.1024 * t2 + 0.0421 * t3
			- 0.0073 * t1 * t2 + 0.024 * t2 * t3 - 0.0118 * t2 * t4
			- 0.0204 * t3 * t4 - 0.008 * t3 * t5 - 0.0241 * t2 * t2
			+ 0.0109 * t4 * t4;
	F[0] = Mass; /* Minimization */
	F[1] = Ain; /* Minimization */
	F[2] = Intrusion; /* Minimization */
	return;
}

/**
 * According to the description in ref (Ganesan et al., 2013)
 * "Swarm intelligence and gravitational search algorithm
 * for multi-objective optimization of synthesis gas production"
 * 3 objectives
 * 3 variables
 */
void Ganesan2013_setup()
{
	gmf_mop.nreal = 3;
	gmf_mop.nobjs = 3;
	gmf_mop.xmin_real = (double*) malloc(sizeof(double) * gmf_mop.nreal);
	gmf_mop.xmax_real = (double*) malloc(sizeof(double) * gmf_mop.nreal);

	gmf_mop.xmin_real[0] = 0.25;
	gmf_mop.xmin_real[1] = 10000.0;
	gmf_mop.xmin_real[2] = 600.0;

	gmf_mop.xmax_real[0] = 0.55;
	gmf_mop.xmax_real[1] = 20000.0;
	gmf_mop.xmax_real[2] = 1100.0;
	return;
}

void Ganesan2013(double *F, double *G, double *xr, int *xi, int *xb)
{
	double O2CH4 = xr[0];
	double GV = xr[1];
	double T = xr[2];

	double HC4_conversion, CO_selectivity, H2_CO_ratio;

	HC4_conversion = (-8.87e-6)
			* (86.74 + 14.6 * O2CH4 - 3.06 * GV + 18.82 * T + 3.14 * O2CH4 * GV
					- 6.91 * O2CH4 * O2CH4 - 13.31 * T * T);

	CO_selectivity = (-2.152e-9)
			* (39.46 + 5.98 * O2CH4 - 2.4 * GV + 13.06 * T + 2.5 * O2CH4 * GV
					+ 1.64 * GV * T - 3.9 * O2CH4 * O2CH4 - 10.15 * T * T
					- 3.69 * GV * GV * O2CH4) + 45.7;

	H2_CO_ratio =
			(4.425e-10)
					* (1.29 - 0.45 * T - 0.112 * O2CH4 * GV - 0.142 * T * GV
							+ 0.109 * O2CH4 * O2CH4 + 0.405 * T * T
							+ 0.167 * T * T * GV) + 0.18;

	F[0] = -HC4_conversion; /* maximization */
	F[1] = -CO_selectivity; /* maximization */
	F[2] = H2_CO_ratio; /* minimization */
	return;
}

/**
 * According to the description in ref (Padhi et al., 2016)
 * "Multi-Objective Optimization of Wire Electrical Discharge Machining (WEDM) Process Parameters
 * Using Weighted Sum Genetic Algorithm Approach"
 * 3 objectives
 * 5 variables
 */
void Padhi2016_setup()
{
	gmf_mop.nobjs = 3;
	gmf_mop.nreal = 5;
	gmf_mop.xmin_real = malloc(sizeof(double) * gmf_mop.nreal);
	gmf_mop.xmax_real = malloc(sizeof(double) * gmf_mop.nreal);

	gmf_mop.xmin_real[0] = 1.0;
	gmf_mop.xmin_real[1] = 10.0;
	gmf_mop.xmin_real[2] = 850.0;
	gmf_mop.xmin_real[3] = 20.0;
	gmf_mop.xmin_real[4] = 4.0;

	gmf_mop.xmax_real[0] = 1.4;
	gmf_mop.xmax_real[1] = 26.0;
	gmf_mop.xmax_real[2] = 1650.0;
	gmf_mop.xmax_real[3] = 40.0;
	gmf_mop.xmax_real[4] = 8.0;
	return;
}
void Padhi2016(double *F, double *G, double *xr, int *xi, int *xb)
{
	double x1 = xr[0];
	double x2 = xr[1];
	double x3 = xr[2];
	double x4 = xr[3];
	double x5 = xr[4];
	double CR;
	double Ra;
	double DD;
	CR = 1.74 + 0.42 * x1 - 0.27 * x2 + 0.087 * x3 - 0.19 * x4 + 0.18 * x5
			+ 0.11 * x1 * x1 + 0.036 * x4 * x4 - 0.025 * x5 * x5
			+ 0.044 * x1 * x2 + 0.034 * x1 * x4 + 0.17 * x1 * x5
			- 0.028 * x2 * x4 + 0.093 * x3 * x4 - 0.033 * x4 * x5;

	Ra = 2.19 + 0.26 * x1 - 0.088 * x2 + 0.037 * x3 - 0.16 * x4 + 0.069 * x5
			+ 0.036 * x1 * x1 + 0.11 * x1 * x3 - 0.077 * x1 * x4
			- 0.075 * x2 * x3 + 0.054 * x2 * x4 + 0.090 * x3 * x5
			+ 0.041 * x4 * x5;

	DD = 0.095 + 0.013 * x1 - 8.625 * 1e-003 * x2 - 5.458 * 1e-003 * x3
			- 0.012 * x4 + 1.462 * 1e-003 * x1 * x1 - 6.635 * 1e-004 * x2 * x2
			- 1.788 * 1e-003 * x4 * x4 - 0.011 * x1 * x2
			- 6.188 * 1e-003 * x1 * x3 + 8.937 * 1e-003 * x1 * x4
			- 4.563 * 1e-003 * x1 * x5 - 0.012 * x2 * x3
			- 1.063 * 1e-003 * x2 * x4 + 2.438 * 1e-003 * x2 * x5
			- 1.937 * 1e-003 * x3 * x4 - 1.188 * 1e-003 * x3 * x5
			- 3.312 * 1e-003 * x4 * x5;

	F[0] = -CR; /* Maximization */
	F[1] = Ra; /* Minimization */
	F[2] = DD; /* Minimization */
	return;
}

/**
 * According to the description in ref (Gao et al., 2020)
 * "Multi-objective optimization of thermal performance of packed bed
 * latent heat thermal storage system based on response surface method"
 * 3 objectives
 * 9 variables
 */
void Gao2020_setup()
{
	gmf_mop.nobjs = 3;
	gmf_mop.nreal = 9;
	gmf_mop.xmin_real = malloc(sizeof(double) * gmf_mop.nreal);
	gmf_mop.xmax_real = malloc(sizeof(double) * gmf_mop.nreal);

	gmf_mop.xmin_real[0] = 40.0;
	gmf_mop.xmin_real[1] = 0.35;
	gmf_mop.xmin_real[2] = 333.0;
	gmf_mop.xmin_real[3] = 20.0;
	gmf_mop.xmin_real[4] = 3000.0;
	gmf_mop.xmin_real[5] = 0.1;
	gmf_mop.xmin_real[6] = 308.0;
	gmf_mop.xmin_real[7] = 150.0;
	gmf_mop.xmin_real[8] = 0.1;

	gmf_mop.xmax_real[0] = 100.0;
	gmf_mop.xmax_real[1] = 0.5;
	gmf_mop.xmax_real[2] = 363.0;
	gmf_mop.xmax_real[3] = 40.0;
	gmf_mop.xmax_real[4] = 4000.0;
	gmf_mop.xmax_real[5] = 3.0;
	gmf_mop.xmax_real[6] = 328.0;
	gmf_mop.xmax_real[7] = 200.0;
	gmf_mop.xmax_real[8] = 2.0;
	return;
}
void Gao2020(double *f, double *g, double *xr, int *xi, int *xb)
{
	double A = xr[0];
	double B = xr[1];
	double C = xr[2];
	double D = xr[3];
	double E = xr[4];
	double F = xr[5];
	double G = xr[6];
	double H = xr[7];
	double J = xr[8];

	double t_eff;
	double Q_eff;
	double Phi_ex;

	t_eff = 171.33 + 23.25 * A - 8.61 * B - 59.85 * C - 66.12 * D - 15.29 * E
			- 83.32 * F + 37.72 * G + 12.67 * H + 0.46 * J - 0.47 * A * B
			- 0.30 * A * C - 6.22 * A * D - 0.62 * A * E - 42.48 * A * F
			+ 3.11 * A * G + 4.45 * A * H - 0.22 * A * J + 7.46 * B * C
			+ 3.28 * B * D + 1.28 * B * E + 1.02 * B * F - 4.02 * B * G
			- 2.29 * B * H - 0.16 * B * J + 19.25 * C * D - 14.83 * C * E
			+ 5.07 * C * F - 37.61 * C * G - 9.11 * C * H - 0.32 * C * J
			+ 8.53 * D * E + 18.46 * D * F - 14.28 * D * G - 7.05 * D * H
			- 0.24 * D * J + 2.05 * E * F + 15.73 * E * G - 0.77 * E * H
			- 0.29 * E * J - 4.77 * F * G + 2.07 * F * H + 0.64 * F * J
			+ 3.41 * G * H + 1.76 * G * J + 0.48 * H * J + 3.64 * A * A
			- 0.99 * B * B + 30.5 * C * C + 21.63 * D * D + 1.72 * E * E
			+ 72.42 * F * F + 11.2 * G * G + 1.86 * H * H - 0.79 * J * J;

	Q_eff = 577.73 - 1.22 * A - 19.56 * B + 102.05 * C - 1.83 * D + 27.28 * E
			+ 2.52 * F + 5.43 * G + 37.48 * H + 0.45 * J + 2.94 * A * B
			- 2.96 * A * C + 0.66 * A * D + 0.09 * A * E - 0.43 * A * F
			+ 0.12 * A * G - 0.43 * A * H - 0.7 * A * J + 8.05 * B * C
			+ 0.53 * B * D + 4.43 * B * E - 0.6 * B * F - 0.46 * B * G
			- 4.97 * B * H + 0.046 * B * J + 0.42 * C * D + 6.03 * C * E
			+ 0.21 * C * F + 2.63 * C * G + 0.17 * C * H - 0.43 * C * J
			+ 6.34 * D * E + 6.36 * D * F + 0.19 * D * G - 0.22 * D * H
			+ 0.39 * D * J - 7.09 * E * F + 3.06 * E * G - 0.15 * E * H
			+ 0.68 * E * J - 0.2 * F * G + 0.14 * F * H + 0.88 * F * J
			+ 0.45 * G * H - 0.014 * G * J + 0.99 * H * J + 0.55 * A * A
			- 4.97 * B * B - 0.47 * C * C - 0.91 * D * D - 2.08 * E * E
			- 1.43 * F * F + 0.43 * G * G + 1.06 * H * H + 0.98 * J * J;

	Phi_ex = 0.81 - 9.26 * 10e-3 * A + 0.014 * B - 0.029 * C - 7.69 * 10e-4 * D
			+ 4.05 * 10e-3 * E + 0.029 * F + 0.075 * G - 0.012 * H
			- 1.04 * 10e-3 * J - 2.63 * 10e-3 * A * B + 1.34 * 10e-4 * A * C
			- 1.48 * 10e-3 * A * D - 7.04 * 10e-4 * A * E + 0.013 * A * F
			+ 6.55 * 10e-4 * A * G - 9.71 * 10e-3 * A * H + 1.08 * 10e-3 * A * J
			+ 2.54 * 10e-3 * B * C - 4.83 * 10e-4 * B * D + 9.63 * 10e-4 * B * E
			+ 1.21 * 10e-3 * B * F - 7.02 * 10e-3 * B * G - 1.21 * 10e-3 * B * H
			+ 1.94 * 10e-5 * B * J - 1.15 * 10e-3 * C * D + 3.60 * 10e-3 * C * E
			+ 5.60 * 10e-3 * C * F - 0.026 * C * G - 4.01 * 10e-3 * C * H
			+ 1.35 * 10e-3 * C * J - 6.93 * 10e-3 * D * E - 3.16 * 10e-3 * D * F
			- 2.38 * 10e-4 * D * G + 7.32 * 10e-4 * D * H + 4.69 * 10e-4 * D * J
			+ 8.18 * 10e-3 * E * F - 5.74 * 10e-3 * E * G + 1.44 * 10e-4 * E * H
			- 9.95 * 10e-5 * E * J - 2.09 * 10e-3 * F * G - 65 * 10e-4 * F * H
			- 1.99 * 10e-3 * F * J + 4.95 * 10e-3 * G * H + 8.70 * 10e-4 * G * J
			+ 4.55 * 10e-4 * H * J - 9.32 * 10e-4 * A * A - 7.61 * 10e-4 * B * B
			+ 0.016 * C * C + 1.24 * 10e-3 * D * D + 9.61 * 10e-4 * E * E
			- 0.024 * F * F - 8.63 * 10e-3 * G * G - 1.90 * 10e-4 * H * H
			- 7.56 * 10e-4 * J * J;

	f[0] = t_eff; /* Minimize */
	f[1] = -Q_eff; /* Maximize */
	f[2] = -Phi_ex; /* Maximize */
	return;
}

/**
 * According to the description in ref (Xu et al., 2020)
 * "Multiobjective Optimization of Milling Parameters for
 * Ultrahigh-Strength Steel AF1410 Based on the NSGA-II Method"
 * 3 objectives
 * 4 variables
 */
void Xu2020_setup()
{
	gmf_mop.nobjs = 3;
	gmf_mop.nreal = 4;
	gmf_mop.xmin_real = malloc(sizeof(double) * gmf_mop.nreal);
	gmf_mop.xmax_real = malloc(sizeof(double) * gmf_mop.nreal);

	gmf_mop.xmin_real[0] = 12.56;
	gmf_mop.xmin_real[1] = 0.02;
	gmf_mop.xmin_real[2] = 1.0;
	gmf_mop.xmin_real[3] = 0.5;

	gmf_mop.xmax_real[0] = 25.12;
	gmf_mop.xmax_real[1] = 0.06;
	gmf_mop.xmax_real[2] = 5.0;
	gmf_mop.xmax_real[3] = 2.0;
	return;
}

void Xu2020(double *F, double *G, double *xr, int *xi, int *xb)
{
	double vc = xr[0];
	double fz = xr[1];
	double ap = xr[2];
	double ae = xr[3];

	double Ft, Ra, MRR;
	double d = 2.5;
	double z = 1.0;

	Ft = -54.3 - 1.18 * vc - 2429 * fz + 104.2 * ap + 129.0 * ae
			- 18.9 * vc * fz - 0.209 * vc * ap - 0.673 * vc * ae + 265 * fz * ap
			+ 1209 * fz * ae + 22.76 * ap * ae + 0.066 * vc * vc
			+ 32117 * fz * fz - 16.98 * ap * ap - 47.6 * ae * ae;
	Ra = 0.227 - 0.0072 * vc + 1.89 * fz - 0.0203 * ap + 0.3075 * ae
			- 0.198 * vc * fz - 0.000955 * vc * ap - 0.00656 * vc * ae
			+ 0.209 * fz * ap + 0.783 * fz * ae + 0.02275 * ap * ae
			+ 0.000355 * vc * vc + 35 * fz * fz + 0.00037 * ap * ap
			- 0.0791 * ae * ae;
	MRR = (1000.0 * vc * fz * z * ap * ae) / (M_PI * d);

	F[0] = Ft; /* Minimization */
	F[1] = Ra; /* Minimization */
	F[2] = -MRR; /* Maximization */
}

/**
 * According to the description in ref (Vaidyanathan et al., 2004)
 * "Computational-fluid-dynamics-based design optimization for single-element
 * rocket injector"
 * 4 objectives
 * 4 variables
 */
void Vaidyanathan2004_setup()
{
	int i;
	gmf_mop.nobjs = 4;
	gmf_mop.nreal = 4;
	gmf_mop.xmin_real = malloc(sizeof(double) * gmf_mop.nreal);
	gmf_mop.xmax_real = malloc(sizeof(double) * gmf_mop.nreal);
	for (i = 0; i < gmf_mop.nreal; ++i)
	{
		gmf_mop.xmin_real[i] = 0.0;
		gmf_mop.xmax_real[i] = 1.0;
	}
	return;
}

void Vaidyanathan2004(double *F, double *G, double *xr, int *xi, int *xb)
{
	double a = xr[0];
	double DHA = xr[1];
	double DOA = xr[2];
	double OPTT = xr[3];

	double Xcc;
	double TFmax;
	double TTmax;
	double TW4;
	// TFmax (minimization)
	TFmax = 0.692 + 0.477 * a - 0.687 * DHA - 0.080 * DOA - 0.0650 * OPTT
			- 0.167 * a * a - 0.0129 * DHA * a + 0.0796 * DHA * DHA
			- 0.0634 * DOA * a - 0.0257 * DOA * DHA + 0.0877 * DOA * DOA
			- 0.0521 * OPTT * a + 0.00156 * OPTT * DHA + 0.00198 * OPTT * DOA
			+ 0.0184 * OPTT * OPTT;

	// Xcc (minimization)
	Xcc = 0.153 - 0.322 * a + 0.396 * DHA + 0.424 * DOA + 0.0226 * OPTT
			+ 0.175 * a * a + 0.0185 * DHA * a - 0.0701 * DHA * DHA
			- 0.251 * DOA * a + 0.179 * DOA * DHA + 0.0150 * DOA * DOA
			+ 0.0134 * OPTT * a + 0.0296 * OPTT * DHA + 0.0752 * OPTT * DOA
			+ 0.0192 * OPTT * OPTT;

	// TW4 (minimization)
	TW4 = 0.758 + 0.358 * a - 0.807 * DHA + 0.0925 * DOA - 0.0468 * OPTT
			- 0.172 * a * a + 0.0106 * DHA * a + 0.0697 * DHA * DHA
			- 0.146 * DOA * a - 0.0416 * DOA * DHA + 0.102 * DOA * DOA
			- 0.0694 * OPTT * a - 0.00503 * OPTT * DHA + 0.0151 * OPTT * DOA
			+ 0.0173 * OPTT * OPTT;

	// TTmax (minimization)
	TTmax = 0.370 - 0.205 * a + 0.0307 * DHA + 0.108 * DOA + 1.019 * OPTT
			- 0.135 * a * a + 0.0141 * DHA * a + 0.0998 * DHA * DHA
			+ 0.208 * DOA * a - 0.0301 * DOA * DHA - 0.226 * DOA * DOA
			+ 0.353 * OPTT * a - 0.0497 * OPTT * DOA - 0.423 * OPTT * OPTT
			+ 0.202 * DHA * a * a - 0.281 * DOA * a * a - 0.342 * DHA * DHA * a
			- 0.245 * DHA * DHA * DOA + 0.281 * DOA * DOA * DHA
			- 0.184 * OPTT * OPTT * a - 0.281 * DHA * a * DOA;

	F[0] = TFmax; /* Minimization */
	F[1] = TW4; /* Minimization */
	F[2] = TTmax; /* Minimization */
	F[3] = Xcc; /* Minimization */
	return;
}

/**
 * According to the description in ref (Chen et al., 2015)
 * "Multiobjective Optimization of Complex Antenna Structures Using Response Surface Models"
 * 5 objectives
 * 6 variables
 */
void Chen2015_setup()
{
	gmf_mop.nobjs = 5;
	gmf_mop.nreal = 6;
	gmf_mop.xmin_real = malloc(sizeof(double) * gmf_mop.nreal);
	gmf_mop.xmax_real = malloc(sizeof(double) * gmf_mop.nreal);
	gmf_mop.xmin_real[0] = 17.5;
	gmf_mop.xmin_real[1] = 17.5;
	gmf_mop.xmin_real[2] = 2.0;
	gmf_mop.xmin_real[3] = 2.0;
	gmf_mop.xmin_real[4] = 5.0;
	gmf_mop.xmin_real[5] = 5.0;

	gmf_mop.xmax_real[0] = 22.5;
	gmf_mop.xmax_real[1] = 22.5;
	gmf_mop.xmax_real[2] = 3.0;
	gmf_mop.xmax_real[3] = 3.0;
	gmf_mop.xmax_real[4] = 7.0;
	gmf_mop.xmax_real[5] = 6.0;
	return;
}
void Chen2015(double *F, double *G, double *xr, int *xi, int *xb)
{
	double l1 = xr[0];
	double w1 = xr[1];
	double l2 = xr[2];
	double w2 = xr[3];
	double a1 = xr[4];
	double b1 = xr[5];
	double a2 = l1 * w1 * l2 * w2;
	double b2 = l1 * w1 * l2 * a1;
	double d2 = w1 * w2 * a1 * b1;

	double F1, F2, F3, F4, F5;

	/* minimization */
	F1 = 502.94 - 27.18 * ((w1 - 20.0) / 0.5) + 43.08 * ((l1 - 20.0) / 2.5)
			+ 47.75 * (a1 - 6.0) + 32.25 * ((b1 - 5.5) / 0.5)
			+ 31.67 * (a2 - 11.0)
			- 36.19 * ((w1 - 20.0) / 0.5) * ((w2 - 2.5) / 0.5)
			- 39.44 * ((w1 - 20.0) / 0.5) * (a1 - 6.0)
			+ 57.45 * (a1 - 6.0) * ((b1 - 5.5) / 0.5);

	F2 = 130.53 + 45.97 * ((l1 - 20.0) / 2.5) - 52.93 * ((w1 - 20.0) / 0.5)
			- 78.93 * (a1 - 6.0) + 79.22 * (a2 - 11.0)
			+ 47.23 * ((w1 - 20.0) / 0.5) * (a1 - 6.0)
			- 40.61 * ((w1 - 20.0) / 0.5) * (a2 - 11.0)
			- 50.62 * (a1 - 6.0) * (a2 - 11.0);

	F3 = 203.16 - 42.75 * ((w1 - 20.0) / 0.5) + 56.67 * (a1 - 6.0)
			+ 19.88 * ((b1 - 5.5) / 0.5) - 12.89 * (a2 - 11.0)
			- 35.09 * (a1 - 6.0) * ((b1 - 5.5) / 0.5)
			- 22.91 * ((b1 - 5.5) / 0.5) * (a2 - 11.0);

	F4 = 0.76 - 0.06 * ((l1 - 20.0) / 2.5) + 0.03 * ((l2 - 2.5) / 0.5)
			+ 0.02 * (a2 - 11.0) - 0.02 * ((b2 - 6.5) / 0.5)
			- 0.03 * ((d2 - 12.0) / 0.5)
			+ 0.03 * ((l1 - 20.0) / 2.5) * ((w1 - 20.0) / 0.5)
			- 0.02 * ((l1 - 20.0) / 2.5) * ((l2 - 2.5) / 0.5)
			+ 0.02 * ((l1 - 20.0) / 2.5) * ((b2 - 6.5) / 0.5);

	/* minimization */
	F5 = 1.08 - 0.12 * ((l1 - 20.0) / 2.5) - 0.26 * ((w1 - 20.0) / 0.5)
			- 0.05 * (a2 - 11.0) - 0.12 * ((b2 - 6.5) / 0.5)
			+ 0.08 * (a1 - 6.0) * ((b2 - 6.5) / 0.5)
			+ 0.07 * (a2 - 6.0) * ((b2 - 5.5) / 0.5);

	F[0] = F1; /* minimization */
	F[1] = -F2; /* maximization */
	F[2] = -F3; /* maximization */
	F[3] = -F4; /* maximization */
	F[4] = F5; /* minimization */
	return;
}

/**
 * According to the description in ref (Ahmad et al., 2017)
 * "Multi-objective optimization in the development of oil and water repellent cellulose
 * fabric based on response surface methodology and the desirability function".
 * NOTE: The sign of the last term is chosen according to the worst scenario
 * for minimization and maximization in the corresponding objective function.
 *
 * 7 objectives
 * 3 variables
 */
void Ahmad2017_setup()
{
	gmf_mop.nobjs = 7;
	gmf_mop.nreal = 3;
	gmf_mop.xmin_real = malloc(sizeof(double) * gmf_mop.nreal);
	gmf_mop.xmax_real = malloc(sizeof(double) * gmf_mop.nreal);
	gmf_mop.xmin_real[0] = 10.0;
	gmf_mop.xmin_real[1] = 10.0;
	gmf_mop.xmin_real[2] = 150.0;

	gmf_mop.xmax_real[0] = 50.0;
	gmf_mop.xmax_real[1] = 50.0;
	gmf_mop.xmax_real[2] = 170.0;
	return;
}
void Ahmad2017(double *F, double *G, double *xr, int *xi, int *xb)
{
	double X1 = xr[0];
	double X2 = xr[1];
	double X3 = xr[2];

	double WCA;
	double OCA;
	double AP;
	double CRA;
	double Stiffness;
	double Tear;
	double Tensile;

	/* maximization */
	WCA = -1331.04 + 1.99 * X1 + 0.33 * X2 + 17.12 * X3 - 0.02 * X1 * X1
			- 0.05 * X3 * X3 - 15.33;

	/* maximization */
	OCA = -4231.14 + 4.27 * X1 + 1.50 * X2 + 52.30 * X3 - 0.04 * X1 * X2
			- 0.04 * X1 * X1 - 0.16 * X3 * X3 - 29.33;

	/* maximization */
	AP = 1766.80 - 32.32 * X1 - 24.56 * X2 - 10.48 * X3 + 0.24 * X1 * X3
			+ 0.19 * X2 * X3 - 0.06 * X1 * X1 - 0.10 * X2 * X2 - 413.33;

	/* maximization */
	CRA = -2342.13 - 1.556 * X1 + 0.77 * X2 + 31.14 * X3 + 0.03 * X1 * X1
			- 0.10 * X3 * X3 - 73.33;

	/* minimization */
	Stiffness = 9.34 + 0.02 * X1 - 0.03 * X2 - 0.03 * X3 - 0.001 * X1 * X2
			+ 0.0009 * X2 * X2 + 0.22;

	/* maximization */
	Tear = 1954.71 + 14.246 * X1 + 5.00 * X2 - 4.30 * X3 - 0.22 * X1 * X1
			- 0.33 * X2 * X2 - 8413.33;

	/* maximization */
	Tensile = 828.16 + 3.55 * X1 + 73.65 * X2 + 10.80 * X3 - 0.56 * X2 * X3
			+ 0.20 * X2 * X2 - 2814.83;

	F[0] = -WCA; /* maximization */
	F[1] = -OCA; /* maximization */
	F[2] = -AP; /* maximization */
	F[3] = -CRA; /* maximization */
	F[4] = Stiffness; /* minimization */
	F[5] = -Tear; /* maximization */
	F[6] = -Tensile; /* maximization */
	return;
}

/** **************************************************************************
 * Setting benchmark
 ** **************************************************************************/
char rwa_name[10][25] =
{ "Subasi2016", "Goel2007", "Liao2008", "Ganesan2013", "Padhi2016", "Gao2020",
		"Xu2020", "Vaidyanathan2004", "Chen2015", "Ahmad2017" };

enum
{
	subasi2016,
	goel2007,
	liao2008,
	ganesan2013,
	padhi2016,
	gao2020,
	xu2020,
	vaidyanathan2004,
	chen2015,
	ahmad2017

};
void (*rwa_mop[10])(double *F, double *G, double *xr, int *xi, int *xb) =
{	Subasi2016,
	Goel2007,
	Liao2008,
	Ganesan2013,
	Padhi2016,
	Gao2020,
	Xu2020,
	Vaidyanathan2004,
	Chen2015,
	Ahmad2017
};

void (*rwa_setup[10])(double *F, double *G, double *xr, int *xi, int *xb) =
{	Subasi2016_setup,
	Goel2007_setup,
	Liao2008_setup,
	Ganesan2013_setup,
	Padhi2016_setup,
	Gao2020_setup,
	Xu2020_setup,
	Vaidyanathan2004_setup,
	Chen2015_setup,
	Ahmad2017_setup
};

void gmf_clean_benchmark()
{
	strcpy(gmf_mop.benchmarck, "");
	strcpy(gmf_mop.name, "");

	gmf_mop.nreal = 0;
	gmf_mop.nint = 0;
	gmf_mop.nbin = 0;
	gmf_mop.ncons = 0;
	gmf_mop.nobjs = 0;
	gmf_mop.dynamic_flag = 0;
	return;
}

/** **************************************************************************
 ** Random numbers
 ** **************************************************************************/
double gmf_rnd_perc(void)
{
	return (double) rand() / (double) RAND_MAX;
}

double rnd_real(double lb, double ub)
{
	assert(lb <= ub);
	return (double) (lb + gmf_rnd_perc() * (ub - lb));
}

void gmf_rwa_rnd_solution(double *xr)
{
	int i;
	for (i = 0; i < gmf_mop.nreal; ++i)
	{
		xr[i] = rnd_real(gmf_mop.xmin_real[i], gmf_mop.xmax_real[i]);
	}
	return;
}

/** **************************************************************************
 ** Settings functions
 ** **************************************************************************/
void gmf_rwa_settings(int function)
{
	void (*gmf_rwa_setup)();
	gmf_clean_benchmark();
	gmf_rwa_setup = rwa_setup[function];
	gmf_rwa_setup();
	gmf_test_problem = rwa_mop[function];

	strcpy(gmf_mop.name, rwa_name[function]);
	strcpy(gmf_mop.benchmarck, "RWA");
	sprintf(gmf_mop.PF_file, " ");
	sprintf(gmf_mop.PS_file, " ");

	gmf_mop.evaluations = 0;
	return;
}

void gmf_rwa_setdown()
{
	free(gmf_mop.xmin_real);
	free(gmf_mop.xmax_real);
	return;
}

/***
 * Display the complete benchmark
 * **/
void gmf_rwa_display_benchmark()
{
	int function;
	printf("RWA:\t");
	for (function = subasi2016; function <= ahmad2017; ++function)
	{
		printf("%s ", rwa_name[function]);
	}
	printf("\n");
	return;
}

int gmf_rwa_setup(char *str_mop)
{
	int function;
	for (function = subasi2016; function <= ahmad2017; ++function)
	{
		if (strcmp(rwa_name[function], str_mop) == 0)
		{
			gmf_rwa_settings(function);
			return 1;
		}
	}
	return 0;
}
