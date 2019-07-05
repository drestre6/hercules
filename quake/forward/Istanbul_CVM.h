/*
 * Istanbul_CVM.h
 *
 *  Created on: Oct 22, 2018
 *      Author: eafit
 */

#ifndef ISTANBUL_CVM_H_
#define ISTANBUL_CVM_H_


/* ===================  pwl_interp_2d_scattered  ===================  */

void material_property_relative_V10_local(double x_input, double y_input, double z_input, double output[3] ) ;

void Istanbul_init ( int32_t myID ) ;
int32_t Istanbul_initparameters ( );


#endif /* ISTANBUL_CVM_H_ */