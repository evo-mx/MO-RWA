/*
 * gmf_dv_rwa.h
 *
 *  Created on: Jun 7, 2022
 *      Author: Saul Zapotecas
 */
#ifndef GMF_DV_RWA_H_
#define GMF_DV_RWA_H_


void gmf_rwa_setdown();
void gmf_rwa_display_benchmark();
int gmf_rwa_setup(char *str_mop);
void gmf_rwa_rnd_solution(double *xr);

#endif /* GMF_DV_RWA_H_ */
