/*
 * Istanbul_CVM.h
 *
 *  Created on: Oct 22, 2018
 *      Author: eafit
 */

#ifndef ISTANBUL_CVM_H_
#define ISTANBUL_CVM_H_


/* ===================  pwl_interp_2d_scattered  ===================  */

double *pwl_interp_2d_scattered_value ( int nd, double xyd[], double zd[], int t_num, int t[], int t_neighbor[], int ni, double xyi[] );

void triangulation_search_delaunay ( int node_num, double node_xy[],
  int triangle_order, int triangle_num, int triangle_node[],
  int triangle_neighbor[], double p[2], int *triangle_index,
  double *alpha, double *beta, double *gamma, int *edge,
  int *step_num );


void material_property_relative_V6(double x_input, double y_input, double z_input, double output[3] );

void material_property_relative_V444S(double x_input, double y_input, double z_input, double output[3] ) ;
void material_property_relative_V444S_nointerp(double x_input, double y_input, double z_input, double output[3] ) ;


#endif /* ISTANBUL_CVM_H_ */
