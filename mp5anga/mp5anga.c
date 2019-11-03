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

//https://farside.ph.utexas.edu/teaching/celestial/Celestialhtml/Celestialhtml.htmlhttps://farside.ph.utexas.edu/teaching/celestial/Celestialhtml/Celestialhtml.html

#include "py/obj.h"
#include "py/runtime.h"
#include "py/builtin.h"
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.1415926535897932
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
#endif
#define floatingpoint double

STATIC floatingpoint __mp5anga_zhr__ = 0.0;
STATIC floatingpoint __mp5anga_latt__ = 0.0;
STATIC floatingpoint __mp5anga_longt__ = 0.0;

floatingpoint rev(floatingpoint x) {
    return  x - floor(x/360.0)*360.0;
}

floatingpoint radians (floatingpoint degrees) {
	return (3.1415926535897932 * degrees / 180.0);
}

floatingpoint degrees (floatingpoint radians) {
	return (180.0 * radians /3.1415926535897932);
}

STATIC mp_obj_t  mp5anga_ts_at_mn(mp_obj_t odd, 
				mp_obj_t omm, 
				mp_obj_t oyy) {
	//https://stjarnhimlen.se/comp/ppcomp.html#3
	//The time scale in these formulae are counted in days. 
	//Hours, minutes, seconds are expressed as fractions of a day. 
	//Day 0.0 occurs at 2000 Jan 0.0 UT (or 1999 Dec 31, 24:00 UT). 
	//This "day number" d is computed as follows (y=year, m=month, D=date, 
	//UT=UT in hours+decimals):
	int dd = mp_obj_get_int (odd);
	int mm = mp_obj_get_int (omm);
	int yy = mp_obj_get_int (oyy);
	int d = 367*yy 
			- 7 * ( yy + (mm+9)/12 ) / 4 
			- 3 * ( ( yy + (mm-9)/7 ) / 100 + 1 ) / 4 
			+ 275*mm/9 + dd - 730515;
	return mp_obj_new_float((floatingpoint)d - (__mp5anga_zhr__)/24.0);//add 2451543.5 to get Julian Date;
}

floatingpoint sun_esrl (floatingpoint* sunriseLSTp, 
			floatingpoint* sunsetLSTp, 
			double julian_day, 
			floatingpoint zhr, 
			floatingpoint latt, 
			floatingpoint longt) {

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
	
	double hour_angle_cosine = cos(radians(90.833))/(cos(radians(latt))*
							cos(radians(sun_declination_deg)))-tan(radians(latt))*
							tan(radians(sun_declination_deg));

	if (hour_angle_cosine > 1.0) {
		//sun never rises
		*sunriseLSTp = 0.0;
		*sunsetLSTp = 0.0;
		return 0.0;
	} else if (hour_angle_cosine < -1.0) {
		//sun never sets
		*sunriseLSTp = 0.0;
		*sunsetLSTp = 0.0;
		return 24.0;
	}

	floatingpoint hour_angle_sunrise_deg = degrees(acos(cos(radians(90.833))/(cos(radians(latt))*
							cos(radians(sun_declination_deg)))-tan(radians(latt))*
							tan(radians(sun_declination_deg)))); //W
	//printf ("hour_angle_sunrise_deg is %f\n", hour_angle_sunrise_deg);

	floatingpoint solar_noon_LST = (720.0-4*longt-equation_of_time_min+zhr*60)/1440.0; //X
	//printf ("solar_noon_LST is %f\n", solar_noon_LST);

	*sunriseLSTp = (solar_noon_LST - (hour_angle_sunrise_deg*4)/1440.0) * 24; //Y

	*sunsetLSTp = (solar_noon_LST + (hour_angle_sunrise_deg*4)/1440.0)*24; //Z

	floatingpoint sunlight_duration_hrs = hour_angle_sunrise_deg * 8 / 60; //AA

	return (sunlight_duration_hrs);
}

floatingpoint nround(double x, unsigned int fac) {
    return round(x*fac)/fac;
}

floatingpoint calculate_moon_sun_long(floatingpoint d, floatingpoint* sLongp) {

	//https://www.stjarnhimlen.se/comp/tutorial.html#7
	//printf ("calculate_moon_sun_long d= %.8f\n", d);

	//------------------------------Sun
	//below are in radians
    floatingpoint w = radians(rev(282.9404 + 4.70935e-5 * d));		//Sun's longitude of perihelion
    floatingpoint Ms = radians(rev(356.0470 + 0.9856002585 * d));	//Sun's mean anomaly
	//floatingpoint oblecl = 23.4393 - 3.563e-7 * d;		//obliquity of the ecliptic

    floatingpoint e = 0.016709 - 1.151e-9 * d;		//eccentricity

	floatingpoint E = Ms + e * sin(Ms) * (1 + e * cos(Ms));//eccentricity anomaly

	//Sun's mean longitude
	floatingpoint Ls = w + Ms;

	//Sun's rectangular coordinates
	floatingpoint x = cos(E) - e;
    floatingpoint y = sin(E) * sqrt(1 - e*e);

	//distance from Sun and true anomaly
	//floatingpoint r = sqrt(x*x + y*y);	//in Earth radii
    floatingpoint v = atan2( y, x );		//true anomaly
	floatingpoint slon = rev(degrees(v + w));
	//printf ("Sun's longitude = %.8f\n", slon);

	//------------------------------Moon
	//all below are in radians
	floatingpoint N = radians(rev(125.1228 - 0.0529538083 * d));   //Longt of ascending node
    const floatingpoint i = 0.089804;			//Inclination in degrees is 5.1454
    w = radians(rev(318.0634 + 0.1643573223 * d));		//Arg. of perigee
    floatingpoint Mm = radians(rev(115.3654 + 13.0649929509 * d));  //Mean eccentricity anomaly

    const floatingpoint a = 60.2666; //Mean distance in Earth equatorial radii
    e = 0.054900;//Eccentricity

	//iterate for accurate eccentricity anomaly
	E = Mm + e * sin(Mm) * (1 + e * cos(Mm));
	floatingpoint eps;
	int iter = 0;
	do {
		eps	= (E - e * sin(E) - Mm) / (1 - e * cos(E));
		E = E - eps;
		if (iter > 50)
			break;
	} while (eps > 1e-5 || eps < -1e-5);

	//compute rectangular (x,y) coordinates in the plane of the lunar orbit
	x = a * (cos(E) - e);
    y = a * sqrt(1 - e*e) * sin(E);

	floatingpoint r = sqrt(x*x + y*y); //distance Earth radii
    v = atan2(y, x); //true anomaly

	floatingpoint xeclip = r * (cos(N) * cos(v+w) - sin(N) * sin(v+w) * cos(i));
    floatingpoint yeclip = r * (sin(N) * cos(v+w) + cos(N) * sin(v+w) * cos(i));
    //floatingpoint zeclip = r * sin(v+w) * sin(i);

	floatingpoint mlon =  rev(degrees(atan2(yeclip, xeclip)));
    //floatingpoint latt  =  atan2(zeclip, sqrt( xeclip*xeclip + yeclip*yeclip));
    //r =  sqrt(xeclip*xeclip + yeclip*yeclip + zeclip*zeclip);

	//Compensate for Moon's perturbations
    //Sun's  mean longitude:        Ls     (already computed as Ls)
    //Moon's mean longitude:        Lm  =  N + w + Mm (for the Moon)
    //Sun's  mean anomaly:          Ms     (already computed as Ms)
    //Moon's mean anomaly:          Mm     (already computed in this function)
    //Moon's mean elongation:       D   =  Lm - Ls
    //Moon's argument of latitude:  F   =  Lm - N
	
	floatingpoint Lm = N + w + Mm;
	floatingpoint D  = Lm - Ls;
	floatingpoint F = Lm - N;
	//printf ("Moon's longt before perturb fix is %f\n", mlon);
	//printf ("Moon's uncorrected ecl. longitude = %.8f\n", mlon);
	mlon += //in degrees
		-1.27388888 * sin(Mm - 2*D)	//Evection -- stjarnhimlen gives -1.274
		+0.65833333 * sin(2*D)		//Variation -- stjarnhimlen give +0.658
		-0.185 * sin(Ms)			//Yearly equation -- stjarnhimlen gives -0.186, but
									//[Chapront-Touzé and Chapront 1988] has 666 arc-seconds
		-0.059 * sin(2*Mm - 2*D)
		-0.057 * sin(Mm - 2*D + Ms)
		+0.053 * sin(Mm + 2*D)
		+0.046 * sin(2*D - Ms)
		+0.041 * sin(Mm - Ms)
		-0.034722222 * sin(D)		//Parallactic equation [Chapront-Touzé and Chapront 1988] has
									//125 arc-seconds = 0.034722222
									//http://www.stjarnhimlen.se/comp/tutorial.html has 0.035
		-0.031 * sin(Mm + Ms)
		-0.015 * sin(2*F - 2*D)		//reduction to the ecliptic from stjarnhimlen -- Wikipedia value is 0.0144
									//stjarnhimlen has 0.015
		+0.011 * sin(Mm - 4*D);
	//printf ("Moon's longt after perturb fix in radians is %f\n", mlon);

	*sLongp = slon;
	//printf ("Sun's ecl. longitude = %.8f\n", slon);
	//printf ("Moon's ecl. longitude = %.8f\n", mlon);
	return mlon;
}

int calculate_tithi_index(floatingpoint d, floatingpoint div) {
	int n;
	floatingpoint mlon, slon;
   	mlon = calculate_moon_sun_long (d, &slon);

	//Calculate Tithi and Paksha
	const floatingpoint fuzz = 0.03;
	n = (int) (nround((rev(mlon-slon+fuzz)/div), 100000));
	//printf ("Diff between Moons and Sun's longitudes = %.8f\n", mlon - slon);
	//printf ("Index of diff between Moons and Sun's longitudes = %.8f and index = %d\n", (mlon - slon)/12.0, n);
	//printf ("Tithi index is n= %d\n", n);
	return n;
}

STATIC mp_obj_t mp5anga_set_zone (mp_obj_t ozhr, mp_obj_t olat, mp_obj_t olong) {
	__mp5anga_latt__ = mp_obj_get_float (olat);
	__mp5anga_longt__ = mp_obj_get_float (olong);
	//ozhr -- time zone w.r.t. UTC
	__mp5anga_zhr__ = mp_obj_get_float(ozhr);
    return mp_const_none;
}

floatingpoint ayanansha(floatingpoint d) {
	floatingpoint t, o, l, ayan;
	
	t = (d+36523.5)/36525;
	o = rev(259.183275-1934.142008333206*t+0.0020777778*t*t);
	l = rev(279.696678+36000.76892*t+0.0003025*t*t);
	ayan = 17.23*sin(radians(o))+1.27*sin(radians(l*2))-(5025.64+1.11*t)*t;
	//Based on Lahiri
	ayan = (ayan-80861.27)/3600.0;

	return ayan;
}

STATIC mp_obj_t mp5anga_tithi (mp_obj_t od) {
	floatingpoint d = mp_obj_get_float (od);
	//printf ("mp5anga_tithi: d is %.8f\n", d);
	int n = calculate_tithi_index (d, 12.0);	
	#ifdef DEBUG
	//printf ("tithi index is %d\n", n);
	#endif
	return mp_obj_new_int(n%30);
}

STATIC mp_obj_t  mp5anga_ayan(mp_obj_t od) {
	floatingpoint ayan = ayanansha(mp_obj_get_float (od));
	return mp_obj_new_float(ayan);
}

STATIC mp_obj_t mp5anga_nakshatra (mp_obj_t od) {
	floatingpoint slon;
	floatingpoint d = mp_obj_get_float (od);
	//return index to nakshatra table
	return mp_obj_new_int((int)(rev(calculate_moon_sun_long (d, &slon) + ayanansha(d))*6/80)%27);
}

STATIC mp_obj_t mp5anga_yoga (mp_obj_t od) {
	floatingpoint mlon, slon;
	floatingpoint d = mp_obj_get_float (od);
	mlon = calculate_moon_sun_long (d, &slon);
	//return index to yoga table
	return mp_obj_new_int((int)(rev(mlon + slon + 2*ayanansha(d))*6/80)%27);
}

STATIC mp_obj_t mp5anga_karana (mp_obj_t od) {
	floatingpoint d = mp_obj_get_float (od);
	int n = calculate_tithi_index (d, 6.0);	
	if(n==0) n=10;
	if(n>=57) n-=50;
	if(n>0 && n<57) n=(n-1)-(floor((n-1)/7)*7); //anirb use floor
	//return index to karana table
	return mp_obj_new_int(n%11);
}

STATIC mp_obj_t mp5anga_rashi (mp_obj_t od) {
	floatingpoint slon;
	floatingpoint d = mp_obj_get_float (od);
	//return index to moon rashi
	return mp_obj_new_int((int)(rev(calculate_moon_sun_long (d, &slon) + ayanansha(d))/30)%12);
}

STATIC mp_obj_t mp5anga_vaara (mp_obj_t od) {
	floatingpoint d = mp_obj_get_float (od) + (__mp5anga_zhr__/24.0);
	//printf ("vaara index: %d\n", (d+5)%7);
	//return index to day
	return mp_obj_new_int(((int)d+5) % 7);
}

floatingpoint fpart (floatingpoint x) {
	double y = x - floor(x);
	if (y < 0)
		y += 1;
	return y;
}

STATIC mp_obj_t mp5anga_sun (mp_obj_t od) {
	char s[26];
	floatingpoint srise, sset, sunlightduration;
	//Bangalore is 12.972N, 77.595E
	floatingpoint d = mp_obj_get_float (od);
	sunlightduration = sun_esrl (&srise, &sset, 
					d,
					__mp5anga_zhr__, 
					__mp5anga_latt__, 
					__mp5anga_longt__);
	if (sunlightduration == 24.0)
		return mp_obj_new_str("0 0 24", 6);
	if (sunlightduration == 0.0)
		return mp_obj_new_str("0 0 0", 5);
	sprintf (s, "%02d:%02d:%02d %02d:%02d:%02d %02d:%02d:%02d", (int)srise, (int)(fpart(srise)*60), (int)nround(fpart(fpart(srise)*60)*60, 1000),
					(int)sset, (int)(fpart(sset)*60), (int)nround(fpart(fpart(sset)*60)*60, 1000),
					(int)sunlightduration, (int)(fpart(sunlightduration)*60), (int)nround(fpart(fpart(sunlightduration)*60)*60, 1000));
	return mp_obj_new_str(s, 26);
}

STATIC MP_DEFINE_CONST_FUN_OBJ_3(mp5anga_set_zone_obj, mp5anga_set_zone);
STATIC MP_DEFINE_CONST_FUN_OBJ_3(mp5anga_ts_at_mn_obj, mp5anga_ts_at_mn);
STATIC MP_DEFINE_CONST_FUN_OBJ_1(mp5anga_tithi_obj, mp5anga_tithi);
STATIC MP_DEFINE_CONST_FUN_OBJ_1(mp5anga_ayan_obj, mp5anga_ayan);
STATIC MP_DEFINE_CONST_FUN_OBJ_1(mp5anga_nakshatra_obj, mp5anga_nakshatra);
STATIC MP_DEFINE_CONST_FUN_OBJ_1(mp5anga_yoga_obj, mp5anga_yoga);
STATIC MP_DEFINE_CONST_FUN_OBJ_1(mp5anga_karana_obj, mp5anga_karana);
STATIC MP_DEFINE_CONST_FUN_OBJ_1(mp5anga_rashi_obj, mp5anga_rashi);
STATIC MP_DEFINE_CONST_FUN_OBJ_1(mp5anga_vaara_obj, mp5anga_vaara);
STATIC MP_DEFINE_CONST_FUN_OBJ_1(mp5anga_sun_obj, mp5anga_sun);

STATIC const mp_rom_map_elem_t mp5anga_module_globals_table[] = {
    { MP_ROM_QSTR(MP_QSTR___name__), MP_ROM_QSTR(MP_QSTR_mp5anga) },
    { MP_ROM_QSTR(MP_QSTR_set_zone), MP_ROM_PTR(&mp5anga_set_zone_obj) },
    { MP_ROM_QSTR(MP_QSTR_ts_at_mn), MP_ROM_PTR(&mp5anga_ts_at_mn_obj) },
    { MP_ROM_QSTR(MP_QSTR_tithi), MP_ROM_PTR(&mp5anga_tithi_obj) },
    { MP_ROM_QSTR(MP_QSTR_ayanansha), MP_ROM_PTR(&mp5anga_ayan_obj) },
    { MP_ROM_QSTR(MP_QSTR_nakshatra), MP_ROM_PTR(&mp5anga_nakshatra_obj) },
    { MP_ROM_QSTR(MP_QSTR_yoga), MP_ROM_PTR(&mp5anga_yoga_obj) },
    { MP_ROM_QSTR(MP_QSTR_karana), MP_ROM_PTR(&mp5anga_karana_obj) },
    { MP_ROM_QSTR(MP_QSTR_rashi), MP_ROM_PTR(&mp5anga_rashi_obj) },
    { MP_ROM_QSTR(MP_QSTR_vaara), MP_ROM_PTR(&mp5anga_vaara_obj) },
    { MP_ROM_QSTR(MP_QSTR_sun), MP_ROM_PTR(&mp5anga_sun_obj) }
};

STATIC MP_DEFINE_CONST_DICT(mp5anga_module_globals, mp5anga_module_globals_table);

const mp_obj_module_t mp5anga_user_cmodule = {
    .base = { &mp_type_module },
    .globals = (mp_obj_dict_t*)&mp5anga_module_globals,
};

MP_REGISTER_MODULE(MP_QSTR_mp5anga, mp5anga_user_cmodule, MODULE_MP5ANGA_ENABLED);


