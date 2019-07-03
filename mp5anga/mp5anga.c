//(c) 2019 Anirban Banerjee
//Licensed under:
//GNU General Public License version 3
/*************************************************************
* Portions of this program are under the following copyright
*************************************************************/
/*
 * Copyright (C) 2012 Santhosh N <santhoshn@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "py/obj.h"
#include "py/runtime.h"
#include "py/builtin.h"
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.141592653
#endif

#ifdef BARE_M
#include "mpconfig.h"
#define sin MICROPY_FLOAT_C_FUN(sin)
#define cos MICROPY_FLOAT_C_FUN(cos)
#define sqrt MICROPY_FLOAT_C_FUN(sqrt)
#define floor MICROPY_FLOAT_C_FUN(floor)
#define atan2 MICROPY_FLOAT_C_FUN(atan2)
#define floatingpoint float
#else
#define floatingpoint double
#endif

#define D2R (M_PI/180.0)
#define R2D (180.0/M_PI)
#define REV(x)	((x)-floor((x)/360.0)*360.0)

STATIC floatingpoint __mp5anga_Ls__, __mp5anga_Lm__, __mp5anga_Ms__, __mp5anga_Mm__ = 0.0;
STATIC floatingpoint __mp5anga_ayan__ = 0.0;
STATIC floatingpoint __mp5anga_mlong__ = 0.0;
STATIC floatingpoint __mp5anga_slong__ = 0.0;
STATIC int __mp5anga_dd__ = 0;
STATIC int __mp5anga_mm__ = 0;
STATIC int __mp5anga_yyyy__ = 0;
STATIC float __mp5anga_zhour__ = 0.0; //local time zone w.r.t. UTC

floatingpoint ayanansha(floatingpoint d) {
	floatingpoint t, o, l, ayan;
	
	t = (d+36523.5)/36525;
	o = 259.183275-1934.142008333206*t+0.0020777778*t*t;
	l = 279.696678+36000.76892*t+0.0003025*t*t;
	ayan = 17.23*sin((o)*D2R)+1.27*sin((l*2)*D2R)-(5025.64+1.11*t)*t;
	//Based on Lahiri
	ayan = (ayan-80861.27)/3600.0;

	return ayan;
}

//Longitude of Sun
floatingpoint slong (floatingpoint d) {
	floatingpoint w, a, e, M, E, x, y, r, v, tmp;

	w = 282.9404+4.70935e-5*d;
	a = 1.000000;
	e = 0.016709-1.151e-9*d;
	M = REV(356.0470+0.9856002585*d);
	__mp5anga_Ms__ = M;
	__mp5anga_Ls__ = w+M;

	tmp = M*D2R;
	E = M+R2D*e*sin(tmp)*(1+e*cos(tmp));

	tmp = E*D2R;
	x = cos(tmp)-e;
	y = sin(tmp)*sqrt(1-e*e);

	r = sqrt(x*x + y*y);
	v = REV(R2D*atan2(y,x));

	return REV(v+w);
}

//Longitude of Moon
floatingpoint mlong(floatingpoint d) {
	floatingpoint N, i, w, a, e, M, E, Et, x, y, r, v, xec, yec, zec, D, F, tmp, tmp1, tmp2, lon; 
	
	N = 125.1228-0.0529538083*d;
 	i = 5.1454;
    	w = REV(318.0634+0.1643573223*d);
    	a = 60.2666;
    	e = 0.054900;
    	M = REV(115.3654+13.0649929509*d);
	__mp5anga_Mm__ = M;
	__mp5anga_Lm__ = N+w+M;

	//Calculate Eccentricity anamoly
	tmp = M*D2R;
	E = M+R2D*e*sin(tmp)*(1+e*cos(tmp));

	tmp = E*D2R;
	Et = E-(E-R2D*e*sin(tmp)-M)/(1-e*cos(tmp));

	do {
		E = Et;
		tmp = E*D2R;
		Et = E-(E-R2D*e*sin(tmp)-M)/(1-e*cos(tmp));
	} while(E-Et>0.005);

	tmp = E*D2R;
	x = a*(cos(tmp)-e);
	y = a*sqrt(1-e*e)*sin(tmp);

	r = sqrt(x*x + y*y);
	v = REV(R2D*atan2(y,x));

	tmp = D2R*N;
	tmp1 = D2R*(v+w);
	tmp2 = D2R*i;
	xec = r*(cos(tmp)*cos(tmp1)-sin(tmp)*sin(tmp1)*cos(tmp2));
	yec = r*(sin(tmp)*cos(tmp1)+cos(tmp)*sin(tmp1)*cos(tmp2));
	zec = r*sin(tmp1)*sin(tmp2);

	//Do some corrections
	D = __mp5anga_Lm__ - __mp5anga_Ls__;
	F = __mp5anga_Lm__ - N;

	lon = R2D*atan2(yec,xec);

	lon+= -1.274*sin((__mp5anga_Mm__-2*D)*D2R);
    lon+= +0.658*sin((2*D)*D2R);
    lon+= -0.186*sin((__mp5anga_Ms__)*D2R);
    lon+= -0.059*sin((2*__mp5anga_Mm__-2*D)*D2R);
    lon+= -0.057*sin((__mp5anga_Mm__-2*D+__mp5anga_Ms__)*D2R);
    lon+= +0.053*sin((__mp5anga_Mm__+2*D)*D2R);
    lon+= +0.046*sin((2*D-__mp5anga_Ms__)*D2R);
    lon+= +0.041*sin((__mp5anga_Mm__-__mp5anga_Ms__)*D2R);
    lon+= -0.035*sin((D)*D2R);
	lon+= -0.031*sin((__mp5anga_Mm__+__mp5anga_Ms__)*D2R);
    lon+= -0.015*sin((2*F-2*D)*D2R);
    lon+= +0.011*sin((__mp5anga_Mm__-4*D)*D2R);

	return REV(lon);
}

STATIC mp_obj_t mp5anga_set_date (mp_obj_t odd, mp_obj_t omm, mp_obj_t oyyyy) {
	__mp5anga_dd__ = mp_obj_get_int(odd);
	__mp5anga_mm__ = mp_obj_get_int(omm);
	__mp5anga_yyyy__ = mp_obj_get_int(oyyyy);
    return mp_const_none;
}

STATIC mp_obj_t mp5anga_set_hour (mp_obj_t ohr, mp_obj_t ozhour) {
	floatingpoint d;
	floatingpoint hr = mp_obj_get_float(ohr);
	__mp5anga_zhour__ = mp_obj_get_float(ozhour);

	//Calculate day number since 2000 Jan 0.0 TDT
	d = (367.0*__mp5anga_yyyy__-7.0*(__mp5anga_yyyy__+(__mp5anga_mm__+9.0)/12.0)/4.0+275.0*__mp5anga_mm__/9.0+__mp5anga_dd__-730530.0);
	__mp5anga_ayan__ = ayanansha (d);
	__mp5anga_slong__ = slong(d + ((hr-__mp5anga_zhour__)/24.0));
	__mp5anga_mlong__ = mlong(d + ((hr-__mp5anga_zhour__)/24.0));
    return mp_const_none;
}

STATIC mp_obj_t mp5anga_tithi (void) {
	floatingpoint tmlon = __mp5anga_mlong__ + ( (__mp5anga_mlong__<__mp5anga_slong__)?360:0 );
	int n = (int)((tmlon - __mp5anga_slong__)/12);
	//return index to tithi and paksha table
	
	#ifdef DEBUG
	printf ("tithi index is %d\n", n);
	#endif
	return mp_obj_new_int(n);
}

STATIC mp_obj_t mp5anga_nakshatra (void) {
	floatingpoint tmlon = REV(__mp5anga_mlong__ + __mp5anga_ayan__);
	int n = (int)(tmlon*6/80);
	#ifdef DEBUG
	printf ("nakshatra index is %d\n", n);
	#endif
	//return index to nakshatra table
	return mp_obj_new_int(n);
}

STATIC mp_obj_t mp5anga_yoga (void) {
	floatingpoint tmlon = __mp5anga_mlong__ + __mp5anga_ayan__;
	floatingpoint tslon = __mp5anga_slong__ + __mp5anga_ayan__;
	#ifdef DEBUG
	printf ("__mp5anga_mlong__ is %f\n", __mp5anga_mlong__);
	printf ("__mp5anga_slong__ is %f\n", __mp5anga_slong__);
	printf ("__mp5anga_ayan__ is %f\n", __mp5anga_ayan__);
	printf ("tmlon is %f\n", tmlon);
	printf ("tslon is %f\n", tslon);
	#endif
	int n = (int)(REV(tmlon+tslon)*6/80);
	#ifdef DEBUG
	printf ("yoga index is %d\n", n);
	#endif
	//return index to yoga table
	return mp_obj_new_int(n);
}

STATIC mp_obj_t mp5anga_karana (void) {
	floatingpoint tmlon = __mp5anga_mlong__ + ((__mp5anga_mlong__ < __mp5anga_slong__)? 360: 0);
	int n = (int)((tmlon - __mp5anga_slong__)/6);
	if (n == 0) n = 10;
	if (n >= 57) n -= 50;
	if (n > 0 && n < 57) n = (n-1) - (floor((n-1)/7)*7); //anirb use floor
	#ifdef DEBUG
	printf ("karana index is %d\n", n);
	#endif

	return mp_obj_new_int(n);
}

STATIC mp_obj_t mp5anga_rashi (void) {
	floatingpoint tmlon = REV(__mp5anga_mlong__ + __mp5anga_ayan__);
	int n = (int)(tmlon/30);
	//return index to moon rashi
	return mp_obj_new_int(n);
}

STATIC mp_obj_t mp5anga_vaara (void) {
	int t1 = (int) floor((__mp5anga_mm__ - 9.0) / 7.0);
	//unsigned int t2 = (int) floor((__mp5anga_yyyy__ + 5001.0 + t1) / 4.0);
	unsigned int t2 = (int) floor((7.0 * (__mp5anga_yyyy__ + 5001.0 + floor((__mp5anga_mm__ - 9.0) / 7.0))) / 4.0);
	unsigned int t3 = (int) floor(275.0 * __mp5anga_mm__ / 9.0);
	unsigned int d = 367.0 * __mp5anga_yyyy__ - t2 + t3 + __mp5anga_dd__ + 1729777.0;
	//printf ("vaara index: %d\n", d%7);
	//return index to day
	return mp_obj_new_int(d % 7);
}

STATIC MP_DEFINE_CONST_FUN_OBJ_3(mp5anga_set_date_obj, mp5anga_set_date);
STATIC MP_DEFINE_CONST_FUN_OBJ_2(mp5anga_set_hour_obj, mp5anga_set_hour);
STATIC MP_DEFINE_CONST_FUN_OBJ_0(mp5anga_tithi_obj, mp5anga_tithi);
STATIC MP_DEFINE_CONST_FUN_OBJ_0(mp5anga_nakshatra_obj, mp5anga_nakshatra);
STATIC MP_DEFINE_CONST_FUN_OBJ_0(mp5anga_yoga_obj, mp5anga_yoga);
STATIC MP_DEFINE_CONST_FUN_OBJ_0(mp5anga_karana_obj, mp5anga_karana);
STATIC MP_DEFINE_CONST_FUN_OBJ_0(mp5anga_rashi_obj, mp5anga_rashi);
STATIC MP_DEFINE_CONST_FUN_OBJ_0(mp5anga_vaara_obj, mp5anga_vaara);

STATIC const mp_rom_map_elem_t mp5anga_module_globals_table[] = {
    { MP_ROM_QSTR(MP_QSTR___name__), MP_ROM_QSTR(MP_QSTR_mp5anga) },
    { MP_ROM_QSTR(MP_QSTR_set_date), MP_ROM_PTR(&mp5anga_set_date_obj) },
    { MP_ROM_QSTR(MP_QSTR_set_hour), MP_ROM_PTR(&mp5anga_set_hour_obj) },
    { MP_ROM_QSTR(MP_QSTR_tithi), MP_ROM_PTR(&mp5anga_tithi_obj) },
    { MP_ROM_QSTR(MP_QSTR_nakshatra), MP_ROM_PTR(&mp5anga_nakshatra_obj) },
    { MP_ROM_QSTR(MP_QSTR_yoga), MP_ROM_PTR(&mp5anga_yoga_obj) },
    { MP_ROM_QSTR(MP_QSTR_karana), MP_ROM_PTR(&mp5anga_karana_obj) },
    { MP_ROM_QSTR(MP_QSTR_rashi), MP_ROM_PTR(&mp5anga_rashi_obj) },
    { MP_ROM_QSTR(MP_QSTR_vaara), MP_ROM_PTR(&mp5anga_vaara_obj) }
};

STATIC MP_DEFINE_CONST_DICT(mp5anga_module_globals, mp5anga_module_globals_table);

const mp_obj_module_t mp5anga_user_cmodule = {
    .base = { &mp_type_module },
    .globals = (mp_obj_dict_t*)&mp5anga_module_globals,
};

MP_REGISTER_MODULE(MP_QSTR_mp5anga, mp5anga_user_cmodule, MODULE_MP5ANGA_ENABLED);


