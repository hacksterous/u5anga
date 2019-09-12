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
#define fmod MICROPY_FLOAT_C_FUN(fmod)
#define asin MICROPY_FLOAT_C_FUN(asin)
#define tan MICROPY_FLOAT_C_FUN(tan)
#define acos MICROPY_FLOAT_C_FUN(acos)
#define floatingpoint float
#else
#define floatingpoint double
#endif

#define D2R (M_PI/180.0)
#define R2D (180.0/M_PI)
#define REV(x)	((x)-floor((x)/360.0)*360.0)

STATIC int __mp5anga_dd__ = 0;
STATIC int __mp5anga_mm__ = 0;
STATIC int __mp5anga_yyyy__ = 0;
STATIC float __mp5anga_hour__; //local time 
STATIC float __mp5anga_zhour__; //time zone w.r.t. UTC

floatingpoint datenum(	int day, 
				int mon, 
				int year, 
				int hour, 
				int min) {
	floatingpoint dNum;
	int cumudays[] = {0, 0,31,59,90,120,151,181,212,243,273,304,334};
	/* Calculate the serial date number:*/
	dNum = (floatingpoint) (365 * year  + cumudays[mon] + day +
		year / 4 - year / 100 + year / 400 +
		(year % 4 != 0) - (year % 100 != 0) + (year % 400 != 0) + 
		((hour * 60 + min) / 1440.0));
	if (mon > 2) {
		if (((year % 4 == 0) && (year % 100 != 0)) || (year % 400 == 0)) {
			dNum += 1.0;
		}
	}
	return (dNum);
}

floatingpoint radians (floatingpoint degrees) {
	return (3.14159265 * degrees / 180.0);
}

floatingpoint degrees (floatingpoint radians) {
	return (180.0 * radians /3.14159265);
}

floatingpoint sun (floatingpoint* sunriseLSTp, 
			floatingpoint* sunsetLSTp, 
			int dd, 
			int mm, 
			int yyyy, 
			floatingpoint zhr, 
			floatingpoint latt, 
			floatingpoint longt) {

	unsigned int day = datenum (dd, mm, yyyy, 12, 0) - 
							datenum (30, 12, 1899, 12, 0);
	//printf ("day is %d\n", day);

	const floatingpoint E = 0.5; //always calculate at noon = 12 hours from local midnight

	floatingpoint julian_day = day + 2415018.5 + E - (zhr/24.0); //F
	//printf ("julian_day is %f\n", julian_day);

	floatingpoint julian_century = (julian_day - 2451545.0)/36525.0; //G
	//printf ("julian_century is %f\n", julian_century);

	floatingpoint sun_geom_mean_long_deg = fmod(280.46646+julian_century*
							(36000.76983 + julian_century*0.0003032), 360.0); //I
	//printf ("sun_geom_mean_long_deg is %f\n", sun_geom_mean_long_deg);

	floatingpoint sun_geom_mean_anom_deg = 357.52911+julian_century*
							(35999.05029 - 0.0001537*julian_century); //J
	//printf ("sun_geom_mean_anom_deg is %f\n", sun_geom_mean_anom_deg);

	floatingpoint earth_orbit_eccentricity = 0.016708634-julian_century*
							(0.000042037+0.0000001267*julian_century); //K
	//printf ("earth_orbit_eccentricity is %f\n", earth_orbit_eccentricity);

	floatingpoint sun_eqn_of_centre = sin(radians(sun_geom_mean_anom_deg))*
							(1.914602-julian_century*(0.004817+0.000014*julian_century))+
							sin(radians(2*sun_geom_mean_anom_deg))*
							(0.019993-0.000101*julian_century)+
							sin(radians(3*sun_geom_mean_anom_deg))*0.000289; //L
	//printf ("sun_eqn_of_centre is %f\n", sun_eqn_of_centre);

	floatingpoint sun_true_long_deg = sun_geom_mean_long_deg + sun_eqn_of_centre; //M
	//printf ("sun_true_long_deg is %f\n", sun_true_long_deg);

	//floatingpoint sun_true_anom_deg = sun_geom_mean_anom_deg + sun_eqn_of_centre; //N
	//printf ("sun_true_anom_deg is %f\n", sun_true_anom_deg);

	//floatingpoint sun_rad_vector_au = (1.000001018*(1-earth_orbit_eccentricity*earth_orbit_eccentricity))/
	//						(1+earth_orbit_eccentricity*cos(radians(sun_true_anom_deg))); //O
	//printf ("sun_rad_vector_au is %f\n", sun_rad_vector_au);

	floatingpoint sun_apparent_long_deg = sun_true_long_deg-0.00569-0.00478*
							sin(radians(125.04-1934.136*julian_century)); //P
	//printf ("sun_apparent_long_deg is %f\n", sun_apparent_long_deg);

	floatingpoint mean_oblique_ecliptic_deg = 23+(26+((21.448-julian_century*(46.815+
							julian_century*(0.00059-julian_century*0.001813))))/60.0)/60.0; //Q
	//printf ("mean_oblique_ecliptic_deg is %f\n", mean_oblique_ecliptic_deg);

	floatingpoint oblique_correction_deg = mean_oblique_ecliptic_deg+
							0.00256*cos(radians(125.04-1934.136*julian_century)); //R
	//printf ("oblique_correction_deg is %f\n", oblique_correction_deg);

	//printf ("radians(sun_apparent_long_deg) is %f\n", radians(sun_apparent_long_deg));

	//floatingpoint sun_right_ascension_deg = degrees(atan2(cos(radians(oblique_correction_deg))*sin(radians(sun_apparent_long_deg)),
	//						cos(radians(sun_apparent_long_deg)))); //S
	//printf ("sun_right_ascension_deg is %f\n", sun_right_ascension_deg);

	floatingpoint sun_declination_deg = degrees(asin(sin(radians(oblique_correction_deg))*
							sin(radians(sun_apparent_long_deg)))); //T
	//printf ("sun_declination_deg is %f\n", sun_declination_deg);

	floatingpoint var_y = tan(radians(oblique_correction_deg/2))*tan(radians(oblique_correction_deg/2)); //U

	floatingpoint equation_of_time_min = 4*degrees(var_y*sin(2*radians(sun_geom_mean_long_deg))-
							2*earth_orbit_eccentricity*sin(radians(sun_geom_mean_anom_deg))+
							4*earth_orbit_eccentricity*var_y*sin(radians(sun_geom_mean_anom_deg))*
							cos(2*radians(sun_geom_mean_long_deg))-
							0.5*var_y*var_y*sin(4*radians(sun_geom_mean_long_deg))-
							1.25*earth_orbit_eccentricity*earth_orbit_eccentricity*
							sin(2*radians(sun_geom_mean_anom_deg))); //V
	//printf ("equation_of_time_min is %f\n", equation_of_time_min);
	
	floatingpoint hour_angle_sunrise_deg = degrees(acos(cos(radians(90.833))/(cos(radians(latt))*
							cos(radians(sun_declination_deg)))-tan(radians(latt))*
							tan(radians(sun_declination_deg)))); //W
	//printf ("hour_angle_sunrise_deg is %f\n", hour_angle_sunrise_deg);

	floatingpoint solar_noon_LST = (720.0-4*longt-equation_of_time_min+zhr*60)/1440.0; //X
	//printf ("solar_noon_LST is %f\n", solar_noon_LST);

	*sunriseLSTp = (solar_noon_LST - (hour_angle_sunrise_deg*4)/1440.0) * 24; //Y

	*sunsetLSTp = (solar_noon_LST + (hour_angle_sunrise_deg*4)/1440.0)*24; //Z

	floatingpoint sunlight_duration_min = hour_angle_sunrise_deg * 8; //AA

	return (sunlight_duration_min);
}

floatingpoint fayanansha(floatingpoint d) {
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
floatingpoint fslong (floatingpoint* Msp, floatingpoint* Lsp, floatingpoint d) {
	floatingpoint w, a, e, M, E, x, y, r, v, tmp;

	w = 282.9404+4.70935e-5*d;
	a = 1.000000;
	e = 0.016709-1.151e-9*d;
	M = REV(356.0470+0.9856002585*d);
	*Msp = M;
	*Lsp = w+M;

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
floatingpoint fmlong(floatingpoint Ms, floatingpoint Ls, floatingpoint d) {
	floatingpoint N, i, w, a, e, M, E, Et, x, y, r, v, xec, yec, zec, D, F, tmp, tmp1, tmp2, lon; 
	floatingpoint Mm, Lm;
	N = 125.1228-0.0529538083*d;
 	i = 5.1454;
    	w = REV(318.0634+0.1643573223*d);
    	a = 60.2666;
    	e = 0.054900;
    	M = REV(115.3654+13.0649929509*d);
	Mm = M;
	Lm = N+w+M;

	//Calculate Eccentricity anomaly
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
	D = Lm - Ls;
	F = Lm - N;

	lon = R2D*atan2(yec,xec);

	lon+= -1.274*sin((Mm-2*D)*D2R);
    lon+= +0.658*sin(2*D*D2R);
    lon+= -0.186*sin(Ms*D2R);
    lon+= -0.059*sin((2*Mm-2*D)*D2R);
    lon+= -0.057*sin((Mm-2*D+Ms)*D2R);
    lon+= +0.053*sin((Mm+2*D)*D2R);
    lon+= +0.046*sin((2*D-Ms)*D2R);
    lon+= +0.041*sin((Mm-Ms)*D2R);
    lon+= -0.035*sin(D*D2R);
	lon+= -0.031*sin((Mm+Ms)*D2R);
    lon+= -0.015*sin((2*F-2*D)*D2R);
    lon+= +0.011*sin((Mm-4*D)*D2R);

	return REV(lon);
}

STATIC mp_obj_t mp5anga_set_date (mp_obj_t odd, mp_obj_t omm, mp_obj_t oyyyy) {
	__mp5anga_dd__ = mp_obj_get_int(odd);
	__mp5anga_mm__ = mp_obj_get_int(omm);
	__mp5anga_yyyy__ = mp_obj_get_int(oyyyy);
    return mp_const_none;
}

STATIC mp_obj_t mp5anga_set_hour (mp_obj_t ohr, mp_obj_t ozhour) {
	if (mp_obj_get_float(ohr) >= 0)
		__mp5anga_hour__ = mp_obj_get_float(ohr);
	if (mp_obj_get_float(ozhour) >= 0)
		__mp5anga_zhour__ = mp_obj_get_float(ozhour);

    return mp_const_none;
}

STATIC mp_obj_t mp5anga_tithi (void) {
	floatingpoint Ms, Ls;
	//Calculate day number since 2000 Jan 0.0 TDT
	floatingpoint d = (367.0*__mp5anga_yyyy__-7.0*(__mp5anga_yyyy__+(__mp5anga_mm__+9.0)/12.0)/4.0+275.0*__mp5anga_mm__/9.0+__mp5anga_dd__-730530.0);
	floatingpoint slong = fslong(&Ms, &Ls, d + ((__mp5anga_hour__ - __mp5anga_zhour__)/24.0));
	floatingpoint mlong = fmlong(Ms, Ls, d + ((__mp5anga_hour__ - __mp5anga_zhour__)/24.0));
	floatingpoint tmlon = mlong + ( (mlong<slong)?360:0 );
	int n = (int)((tmlon - slong)/12);
	//return index to tithi and paksha table
	
	#ifdef DEBUG
	printf ("tithi index is %d\n", n);
	#endif
	return mp_obj_new_int(n);
}

STATIC mp_obj_t mp5anga_nakshatra (void) {
	floatingpoint Ms, Ls;
	//Calculate day number since 2000 Jan 0.0 TDT
	floatingpoint d = (367.0*__mp5anga_yyyy__-7.0*(__mp5anga_yyyy__+(__mp5anga_mm__+9.0)/12.0)/4.0+275.0*__mp5anga_mm__/9.0+__mp5anga_dd__-730530.0);
	floatingpoint ayanansha = fayanansha (d);
	fslong(&Ms, &Ls, d + ((__mp5anga_hour__ - __mp5anga_zhour__)/24.0));
	floatingpoint mlong = fmlong(Ms, Ls, d + ((__mp5anga_hour__ - __mp5anga_zhour__)/24.0));
	floatingpoint tmlon = REV(mlong + ayanansha);
	int n = (int)(tmlon*6/80);
	#ifdef DEBUG
	printf ("mp5anga_nakshatra: d is %f\n", d);
	printf ("mp5anga_nakshatra: ayanansha is %f\n", ayanansha);
	printf ("mp5anga_nakshatra: mlong is %f\n", mlong);
	printf ("mp5anga_nakshatra: Ms is %f\n", Ms);
	printf ("mp5anga_nakshatra: Ls is %f\n", Ls);
	printf ("nakshatra index is %d\n", n);
	#endif
	//return index to nakshatra table
	return mp_obj_new_int(n);
}

STATIC mp_obj_t mp5anga_yoga (void) {
	floatingpoint Ms, Ls;
	//Calculate day number since 2000 Jan 0.0 TDT
	floatingpoint d = (367.0*__mp5anga_yyyy__-7.0*(__mp5anga_yyyy__+(__mp5anga_mm__+9.0)/12.0)/4.0+275.0*__mp5anga_mm__/9.0+__mp5anga_dd__-730530.0);
	floatingpoint ayanansha = fayanansha (d);
	floatingpoint tslon = fslong(&Ms, &Ls, d + ((__mp5anga_hour__ - __mp5anga_zhour__)/24.0)) + ayanansha;
	floatingpoint tmlon = fmlong(Ms, Ls, d + ((__mp5anga_hour__ - __mp5anga_zhour__)/24.0)) + ayanansha;
	int n = (int)(REV(tmlon+tslon)*6/80);
	#ifdef DEBUG
	printf ("tslon is %f\n", tslon);
	printf ("tmlontmlon is %f\n", tmlon);
	printf ("ayanansha is %f\n", ayanansha);
	printf ("yoga index is %d\n", n);
	#endif
	//return index to yoga table
	return mp_obj_new_int(n);
}

STATIC mp_obj_t mp5anga_karana (void) {
	floatingpoint Ms, Ls;
	//Calculate day number since 2000 Jan 0.0 TDT
	floatingpoint d = (367.0*__mp5anga_yyyy__-7.0*(__mp5anga_yyyy__+(__mp5anga_mm__+9.0)/12.0)/4.0+275.0*__mp5anga_mm__/9.0+__mp5anga_dd__-730530.0);
	floatingpoint slong = fslong(&Ms, &Ls, d + ((__mp5anga_hour__ - __mp5anga_zhour__)/24.0));
	floatingpoint mlong = fmlong(Ms, Ls, d + ((__mp5anga_hour__ - __mp5anga_zhour__)/24.0));

	floatingpoint tmlon = mlong + ((mlong < slong)? 360: 0);
	int n = (int)((tmlon - slong)/6);
	if (n == 0) n = 10;
	if (n >= 57) n -= 50;
	if (n > 0 && n < 57) n = (n-1) - (floor((n-1)/7)*7); //anirb use floor
	#ifdef DEBUG
	printf ("karana index is %d\n", n);
	#endif

	return mp_obj_new_int(n);
}

STATIC mp_obj_t mp5anga_rashi (void) {
	floatingpoint Ms, Ls;
	//Calculate day number since 2000 Jan 0.0 TDT
	floatingpoint d = (367.0*__mp5anga_yyyy__-7.0*(__mp5anga_yyyy__+(__mp5anga_mm__+9.0)/12.0)/4.0+275.0*__mp5anga_mm__/9.0+__mp5anga_dd__-730530.0);
	floatingpoint ayanansha = fayanansha (d);
	fslong(&Ms, &Ls, d + ((__mp5anga_hour__ - __mp5anga_zhour__)/24.0));
	floatingpoint mlong = fmlong(Ms, Ls, d + ((__mp5anga_hour__ - __mp5anga_zhour__)/24.0));

	floatingpoint tmlon = REV(mlong + ayanansha);
	int n = (int)(tmlon/30);
	//return index to moon rashi
	return mp_obj_new_int(n);
}

STATIC mp_obj_t mp5anga_vaara (void) {
	//unsigned int t1 = (int) floor((7.0 * (__mp5anga_yyyy__ + 5001.0 + floor((__mp5anga_mm__ - 9.0) / 7.0))) / 4.0);
	//unsigned int t2 = (int) floor(275.0 * __mp5anga_mm__ / 9.0);
	//unsigned int d = 367.0 * __mp5anga_yyyy__ - t1 + t2 + __mp5anga_dd__ + 1729777.0; // --- bug!
	unsigned int d = datenum (__mp5anga_dd__, __mp5anga_mm__, __mp5anga_yyyy__, 12, 0);
	//printf ("vaara index: %d\n", (d+5)%7);
	//return index to day
	return mp_obj_new_int((d+5) % 7);
}

STATIC mp_obj_t mp5anga_srise (mp_obj_t olat, mp_obj_t olong) {
	floatingpoint srise, sset;
	//Bangalore is 12.972N, 77.595E
	floatingpoint latitude = mp_obj_get_float (olat);
	floatingpoint longitude = mp_obj_get_float (olong);
	sun (&srise, &sset, 
					__mp5anga_dd__, __mp5anga_mm__, __mp5anga_yyyy__, 
					__mp5anga_zhour__, latitude, 
					longitude);
	return mp_obj_new_float(srise);
}

STATIC mp_obj_t mp5anga_sset (mp_obj_t olat, mp_obj_t olong) {
	floatingpoint srise, sset;
	floatingpoint latitude = mp_obj_get_float (olat);
	floatingpoint longitude = mp_obj_get_float (olong);
	sun (&srise, &sset, 
					__mp5anga_dd__, __mp5anga_mm__, __mp5anga_yyyy__, 
					__mp5anga_zhour__, latitude, 
					longitude);

	return mp_obj_new_float(sset);
}

STATIC mp_obj_t mp5anga_shrs (mp_obj_t olat, mp_obj_t olong) {
	floatingpoint srise, sset;
	floatingpoint latitude = mp_obj_get_float (olat);
	floatingpoint longitude = mp_obj_get_float (olong);
	floatingpoint sunhrs = sun (&srise, &sset, 
					__mp5anga_dd__, __mp5anga_mm__, __mp5anga_yyyy__, 
					__mp5anga_zhour__, latitude, 
					longitude);

	return mp_obj_new_float(sunhrs);
}

STATIC MP_DEFINE_CONST_FUN_OBJ_3(mp5anga_set_date_obj, mp5anga_set_date);
STATIC MP_DEFINE_CONST_FUN_OBJ_2(mp5anga_set_hour_obj, mp5anga_set_hour);
STATIC MP_DEFINE_CONST_FUN_OBJ_0(mp5anga_tithi_obj, mp5anga_tithi);
STATIC MP_DEFINE_CONST_FUN_OBJ_0(mp5anga_nakshatra_obj, mp5anga_nakshatra);
STATIC MP_DEFINE_CONST_FUN_OBJ_0(mp5anga_yoga_obj, mp5anga_yoga);
STATIC MP_DEFINE_CONST_FUN_OBJ_0(mp5anga_karana_obj, mp5anga_karana);
STATIC MP_DEFINE_CONST_FUN_OBJ_0(mp5anga_rashi_obj, mp5anga_rashi);
STATIC MP_DEFINE_CONST_FUN_OBJ_0(mp5anga_vaara_obj, mp5anga_vaara);
STATIC MP_DEFINE_CONST_FUN_OBJ_2(mp5anga_srise_obj, mp5anga_srise);
STATIC MP_DEFINE_CONST_FUN_OBJ_2(mp5anga_sset_obj, mp5anga_sset);
STATIC MP_DEFINE_CONST_FUN_OBJ_2(mp5anga_shrs_obj, mp5anga_shrs);

STATIC const mp_rom_map_elem_t mp5anga_module_globals_table[] = {
    { MP_ROM_QSTR(MP_QSTR___name__), MP_ROM_QSTR(MP_QSTR_mp5anga) },
    { MP_ROM_QSTR(MP_QSTR_set_date), MP_ROM_PTR(&mp5anga_set_date_obj) },
    { MP_ROM_QSTR(MP_QSTR_set_hour), MP_ROM_PTR(&mp5anga_set_hour_obj) },
    { MP_ROM_QSTR(MP_QSTR_tithi), MP_ROM_PTR(&mp5anga_tithi_obj) },
    { MP_ROM_QSTR(MP_QSTR_nakshatra), MP_ROM_PTR(&mp5anga_nakshatra_obj) },
    { MP_ROM_QSTR(MP_QSTR_yoga), MP_ROM_PTR(&mp5anga_yoga_obj) },
    { MP_ROM_QSTR(MP_QSTR_karana), MP_ROM_PTR(&mp5anga_karana_obj) },
    { MP_ROM_QSTR(MP_QSTR_rashi), MP_ROM_PTR(&mp5anga_rashi_obj) },
    { MP_ROM_QSTR(MP_QSTR_vaara), MP_ROM_PTR(&mp5anga_vaara_obj) },
    { MP_ROM_QSTR(MP_QSTR_srise), MP_ROM_PTR(&mp5anga_srise_obj) },
    { MP_ROM_QSTR(MP_QSTR_sset), MP_ROM_PTR(&mp5anga_sset_obj) },
    { MP_ROM_QSTR(MP_QSTR_shrs), MP_ROM_PTR(&mp5anga_shrs_obj) }
};

STATIC MP_DEFINE_CONST_DICT(mp5anga_module_globals, mp5anga_module_globals_table);

const mp_obj_module_t mp5anga_user_cmodule = {
    .base = { &mp_type_module },
    .globals = (mp_obj_dict_t*)&mp5anga_module_globals,
};

MP_REGISTER_MODULE(MP_QSTR_mp5anga, mp5anga_user_cmodule, MODULE_MP5ANGA_ENABLED);


