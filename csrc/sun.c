/*
 * Copyright (C) 2019 Anirban Banerjee <anirbax@gmail.com>
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

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>

typedef struct {
	char today[16];
	char vaara[16];
	char yoga[64];
	char nakshatra[64];
	char tithi[64];
	char karana[64];
	char paksha[64];
	char rashi[32];
} panchanga;

//static char __g_month[][15]={"January","February","March","April","May","June",
//		  "July","August","September","October","November","December"};
//
//static char __g_rashi[][15]={"Mesha","Vrishabha","Mithuna","Karka","Simha","Kanya","Tula",
//		   "Vrischika","Dhanu","Makara","Kumbha","Meena"};

static char __g_day[][15]={"Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday"};

static char __g_tithi[][15]={"Prathamaa","Dvitiya","Tritiya","Chaturthi","Panchami",
		  "Shashthi","Saptami","Ashtami","Navami","Dashami","Ekadashi",
		  "Dvadashi","Trayodashi","Chaturdashi","Purnima","Pratipada",
		  "Dvitiya","Tritiya","Chaturthi","Panchami","Shashthi",
		  "Saptami","Ashtami","Navami","Dashami","Ekadashi","Dvadashi",
		  "Trayodashi","Chaturdashi","Amaavasya"};

//static char __g_karan[][15]={"Bava","Baalava","Kaulava","Taitula","Garija","Vanija",
//		  "Vishti","Shakuni","Chatushpada","Naga","Kimstughna"};
//
//static char __g_yoga[][15]={"Vishakumbha","Priti","Ayushman","Saubhagya","Shobhana",
//		 "Atiganda","Sukarman","Dhriti","Shula","Ganda","Vriddhi",
//		 "Dhruva","Vyaghata","Harshana","Vajra","Siddhi","Vyatipata",
//		 "Variyan","Parigha","Shiva","Siddha","Saadhya","Shubha","Shukla",
//		 "Brahma","Indra","Vaidhriti"};
//
//static char __g_nakshatra[][20]={"Ashvini","Bharani","Krittika","Rohini","Mrigashira","Ardra",
//		      "Punarvasu","Pushya","Ashlesa","Magha","Purva Phalguni","Uttara Phalguni",
//		      "Hasta","Chitra","Svaati","Vishakha","Anuradha","Jyeshtha","Mula",
//		      "Purva Ashadha","Uttara Ashadha","Shravana","Dhanishtha","Shatabhisha",
//		      "Purva Bhaadra","Uttara Bhaadra","Revati"};

#define PI          3.14159265358979323846
#define RADEG       (180.0/PI)
#define DEGRAD      (PI/180.0)
#define sind(x)     sin((x)*DEGRAD)
#define cosd(x)     cos((x)*DEGRAD)
#define tand(x)     tan((x)*DEGRAD)
#define asind(x)    (RADEG*asin(x))
#define acosd(x)    (RADEG*acos(x))
#define atand(x)    (RADEG*atan(x))
#define atan2d(y,x) (RADEG*atan2((y),(x)))


double rev(double x) {
    return  x - floor(x/360.0)*360.0;
}

double datenum(	int day, 
				int mon, 
				int year, 
				int hour, 
				int min, 
				int sec) {
	int tmp1, tmp2, tmp3;
	double	tmp4, tmp5;
	double dNum;
	int cumudays[] = {0, 0,31,59,90,120,151,181,212,243,273,304,334};
	/* Calculate the serial date number:*/
	tmp1 = 365 * year  + cumudays[mon] + day;
	tmp2 = year / 4 - year / 100 + year / 400;
	tmp3 = (year % 4 != 0) - (year % 100 != 0) + (year % 400 != 0);
	tmp4 = (double) (tmp1+tmp2+tmp3);
	tmp5 = (hour * 3600 + min * 60 + sec) / 86400.0;
	dNum = tmp4 + tmp5;
	if (mon > 2) {
		if (((year % 4 == 0) && (year % 100 != 0)) || (year % 400 == 0)) {
			dNum += 1.0;
		}
	}
	return (dNum);
}

double timescale(int dd, 
				int mm, 
				int yy, 
				float hr, 
				float zhr) {
	//https://stjarnhimlen.se/comp/ppcomp.html#3
	//The time scale in these formulae are counted in days. 
	//Hours, minutes, seconds are expressed as fractions of a day. 
	//Day 0.0 occurs at 2000 Jan 0.0 UT (or 1999 Dec 31, 24:00 UT). 
	//This "day number" d is computed as follows (y=year, m=month, D=date, 
	//UT=UT in hours+decimals):
	int d = 367*yy 
			- 7 * ( yy + (mm+9)/12 ) / 4 
			- 3 * ( ( yy + (mm-9)/7 ) / 100 + 1 ) / 4 
			+ 275*mm/9 + dd - 730515;

	return (double) d + (hr - zhr)/24.0;//add 2451543.5 to get Julian Date;
}

double radians (double degrees) {
	return (3.14159265358979323846 * degrees / 180.0);
}

double degrees (double radians) {
	return (180.0 * radians /3.14159265358979323846);
}

int sign (double y) {
	return ((int) (y > 0) - (int) (y < 0));
}

double ipart (double x) {
	double ix = floor(fabs(x));
	if (ix < 0)
		ix += 1;
	ix = ix * (double) sign(x);
	return (ix);
}

double fpart (double x) {
	double y = x - floor(x);
	if (y < 0)
		y += 1;
	return y;
}

double JulianDay (int date, int month, int year, double hr, double zhr) {
	//http://www.geoastro.de/elevazmoon/basics/meeus.htm -- this has a bug
	//----------------valid from 1900/3/1 to 2100/2/28
    if (month<=2) {month=month+12; year=year-1;}

	//Don't use the code given in http://www.geoastro.de/elevazmoon/basics/meeus.htm
    //return (365.25*year) + floor(30.6001*(month+1)) - 15.0 + 1720996.5 + date + UT/24.0;
	
	//Use formula given in 
	//A Physically-Based Night Sky ModelA Physically-Based Night Sky Model
	//by Julie Dorsey
	/*
	JD = 1720996.5 − floor(year/100) + floor(year/400) + 
			floor(365.25 * year) + floor (30.6001 * (month + 1)) + 
			D + (h + (m + s/60)/60)/24
			### h:m:s is local time with zone correction
	NOTE: deltaT of 65 seconds (=0.018055555555556hr) has not been added to hr
	*/
    return floor(year/400) - floor (year/100) + floor(365.25*year) + 
			floor(30.6001*(month+1)) + 1720996.5 + date + (hr - zhr)/24.0;
}

double JulianDay2 (int dd, 
			int mm, 
			int yyyy, 
			double hr,
			double zhr) {

	//NOTE: deltaT of 65 seconds (=0.018055555555556hr) has not been added to hr
	double day = datenum (dd, mm, yyyy, 0, 0, 0) - 
							datenum (30, 12, 1899, 0, 0, 0);
	//printf ("A: day is %.15f\n", day);

	//const double E = 0.5; //always calculate at noon = 12 hours from local midnight
	double E = hr/24.0;

	double julian_day = day + 2415018.5 + E - (zhr/24.0);
	//printf ("B: JulianDay2 is %.15f\n", julian_day);

	return julian_day;
	
}

//Modified Julian Day
double ModifiedJulianDay (int dd, int mm, int yy, double h, double z) {
	//http://www.stargazing.net/kepler/riset.bas
	//returns modified julian date = Julian Day - 2400000.5
	//number of days since 1858 Nov 17 00:00h
	//valid for any date since 4713 BC
	//assumes gregorian calendar after 1582 Oct 15, Julian before
	//Years BC assumed in calendar format, i.e. the year before 1 AD is 1 BC
	int d = dd;
	int m = mm;
	int y = yy;
	double b;
	double a = 10000.0 * y + 100.0 * m + d;
	//printf ("------ModifiedJulianDay: b is %.15f\n", b);

	if (y < 0)
		y = y + 1;
	if (m <= 2) {
		m = m + 12;
		y = y - 1;
	}
	if (a <= 15821004.1)
	   b = -2 + ipart((y + 4716) / 4) - 1179;
	else
	   b = ipart(y / 400) - ipart(y / 100) + ipart(y / 4);
	//printf ("------ModifiedJulianDay: b is %.15f\n", b);
	a = 365.0 * y - 679004.0;
	a = a + b + ipart(30.6001 * (m + 1)) + d + (h/24) - (z/24);
	//NOTE: deltaT of 65 seconds (=0.018055555555556hr) has not been added to hr
	//printf("Modified Julian Day = %.15f\n", a);
	//printf("Julian Day 3 = %.15f\n", a+2400000.5);
	return (a);
}

double sun_esrl (double* sunriseLSTp, 
			double* sunsetLSTp, 
			int dd, 
			int mm, 
			int yyyy, 
			double zhr, 
			double latt, 
			double longt, int altfn) {

	// Reverse engineered from the NOAA Excel:
	//https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html

	double julian_day;
	if (altfn)
		julian_day = JulianDay (dd, mm, yyyy, 12.0, zhr);
	else
		julian_day = JulianDay2 (dd, mm, yyyy, 12.0, zhr);
	printf ("sun_esrl: julian_day is %.15f\n", julian_day);

	double julian_century = (julian_day - 2451545.0)/36525.0; //G
	//printf ("julian_century is %.15f\n", julian_century);

	double sun_geom_mean_long_deg = fmod(280.46646+julian_century*
							(36000.76983 + julian_century*0.0003032), 360.0); //I
	//printf ("sun_geom_mean_long_deg is %.15f\n", sun_geom_mean_long_deg);

	double sun_geom_mean_anom_deg = 357.52911+julian_century*
							(35999.05029 - 0.0001537*julian_century); //J
	//printf ("sun_geom_mean_anom_deg is %.15f\n", sun_geom_mean_anom_deg);

	double earth_orbit_eccentricity = 0.016708634-julian_century*
							(0.000042037+0.0000001267*julian_century); //K
	//printf ("earth_orbit_eccentricity is %.15f\n", earth_orbit_eccentricity);

	double sun_eqn_of_centre = sin(radians(sun_geom_mean_anom_deg))*
							(1.914602-julian_century*(0.004817+0.000014*julian_century))+
							sin(radians(2*sun_geom_mean_anom_deg))*
							(0.019993-0.000101*julian_century)+
							sin(radians(3*sun_geom_mean_anom_deg))*0.000289; //L
	//printf ("sun_eqn_of_centre is %.15f\n", sun_eqn_of_centre);

	double sun_true_long_deg = sun_geom_mean_long_deg + sun_eqn_of_centre; //M
	//printf ("sun_true_long_deg is %.15f\n", sun_true_long_deg);

	double sun_true_anom_deg = sun_geom_mean_anom_deg + sun_eqn_of_centre; //N
	//printf ("sun_true_anom_deg is %.15f\n", sun_true_anom_deg);

	double sun_rad_vector_au = (1.000001018*(1-earth_orbit_eccentricity*earth_orbit_eccentricity))/
							(1+earth_orbit_eccentricity*cos(radians(sun_true_anom_deg))); //O
	//printf ("sun_rad_vector_au is %.15f\n", sun_rad_vector_au);

	double sun_apparent_long_deg = sun_true_long_deg-0.00569-0.00478*
							sin(radians(125.04-1934.136*julian_century)); //P
	//printf ("sun_apparent_long_deg is %.15f\n", sun_apparent_long_deg);

	double mean_oblique_ecliptic_deg = 23+(26+((21.448-julian_century*(46.815+
							julian_century*(0.00059-julian_century*0.001813))))/60.0)/60.0; //Q
	//printf ("mean_oblique_ecliptic_deg is %.15f\n", mean_oblique_ecliptic_deg);

	double oblique_correction_deg = mean_oblique_ecliptic_deg+
							0.00256*cos(radians(125.04-1934.136*julian_century)); //R
	//printf ("oblique_correction_deg is %.15f\n", oblique_correction_deg);

	//printf ("radians(sun_apparent_long_deg) is %.15f\n", radians(sun_apparent_long_deg));

	double sun_right_ascension_deg = degrees(atan2(cos(radians(oblique_correction_deg))*sin(radians(sun_apparent_long_deg)),
							cos(radians(sun_apparent_long_deg)))); //S
	//printf ("sun_right_ascension_deg is %.15f\n", sun_right_ascension_deg);

	double sun_declination_deg = degrees(asin(sin(radians(oblique_correction_deg))*
							sin(radians(sun_apparent_long_deg)))); //T
	//printf ("sun_declination_deg is %.15f\n", sun_declination_deg);

	double var_y = tan(radians(oblique_correction_deg/2))*tan(radians(oblique_correction_deg/2)); //U

	double equation_of_time_min = 4*degrees(var_y*sin(2*radians(sun_geom_mean_long_deg))-
							2*earth_orbit_eccentricity*sin(radians(sun_geom_mean_anom_deg))+
							4*earth_orbit_eccentricity*var_y*sin(radians(sun_geom_mean_anom_deg))*
							cos(2*radians(sun_geom_mean_long_deg))-
							0.5*var_y*var_y*sin(4*radians(sun_geom_mean_long_deg))-
							1.25*earth_orbit_eccentricity*earth_orbit_eccentricity*
							sin(2*radians(sun_geom_mean_anom_deg))); //V
	//printf ("equation_of_time_min is %.15f\n", equation_of_time_min);
	
	double hour_angle_sunrise_deg = degrees(acos(cos(radians(90.833))/(cos(radians(latt))*
							cos(radians(sun_declination_deg)))-tan(radians(latt))*
							tan(radians(sun_declination_deg)))); //W
	//printf ("hour_angle_sunrise_deg is %.15f\n", hour_angle_sunrise_deg);

	double solar_noon_LST = (720.0-4*longt-equation_of_time_min+zhr*60)/1440.0; //X
	//printf ("solar_noon_LST is %.15f\n", solar_noon_LST);

	*sunriseLSTp = (solar_noon_LST - (hour_angle_sunrise_deg*4)/1440.0) * 24; //Y

	*sunsetLSTp = (solar_noon_LST + (hour_angle_sunrise_deg*4)/1440.0)*24; //Z

	double sunlight_duration_min = hour_angle_sunrise_deg * 8; //AA

	return (sunlight_duration_min);
}

void calculate_tithi(int dd, int mm, int yy, double hr, double zhr, panchanga *pdata) {
	//https://www.stjarnhimlen.se/comp/tutorial.html#7
	double d = timescale (dd, mm, yy, hr, zhr);

	//------------------------------Sun
	//below are in degrees
    double w = rev(282.9404 + 4.70935e-5 * d);		//Sun's longitude of perihelion
    double Ms = rev(356.0470 + 0.9856002585 * d);	//Sun's mean anomaly
	//double oblecl = 23.4393 - 3.563e-7 * d;		//obliquity of the ecliptic

    double e = 0.016709 - 1.151e-9 * d;		//eccentricity

	double E = Ms + (180/PI) * e * sind(Ms) * (1 + e * cosd(Ms));//eccentricity anomaly

	//Sun's mean longitude
	double Ls = rev(w + Ms);
	//printf ("Sun's eccentricity anomaly E = %.8f\n", E);

	//Sun's rectangular coordinates
	double x = cosd(E) - e;
    double y = sind(E) * sqrt(1 - e*e);

	//distance from Sun and true anomaly
	//double r = sqrt(x*x + y*y);	//in Earth radii
    double v = atan2d( y, x );	//true anomaly
	double slon = rev(v + w);
	printf ("Sun's longitude = %.8f\n", slon);

	//------------------------------Moon
	//all below are in degrees
	double N = rev(125.1228 - 0.0529538083  * d);   //Longt of ascending node
    const double i = 5.1454;						//Inclination
    w = rev(318.0634 + 0.1643573223  * d);   //Arg. of perigee
    double Mm = rev(115.3654 + 13.0649929509 * d);  //Mean eccentricity anomaly

    const double a = 60.2666; //Mean distance in Earth equatorial radii
    e = 0.054900;//Eccentricity

	//iterate for accurate eccentricity anomaly
	E = Mm + (180/PI) * e * sind(Mm) * (1 + e * cosd(Mm));
	double eps;
	do {
		eps	= (E - (180/PI) * e * sind(E) - Mm) / (1 - e * cosd(E));
		E = E - eps;
	} while (eps > 0.005 || eps < -0.005);
	//compute rectangular (x,y) coordinates in the plane of the lunar orbit
	x = a * (cosd(E) - e);
    y = a * sqrt(1 - e*e) * sind(E);
	//printf ("Moon x = %.8f\n", x);
	//printf ("Moon y = %.8f\n", y);

	double r = sqrt(x*x + y*y); //distance Earth radii
    v = atan2d(y, x); //true anomaly

	double xeclip = r * (cosd(N) * cosd(v+w) - sind(N) * sind(v+w) * cosd(i));
    double yeclip = r * (sind(N) * cosd(v+w) + cosd(N) * sind(v+w) * cosd(i));
    //double zeclip = r * sind(v+w) * sind(i);
	//printf ("Moon xeclip = %.8f\n", xeclip);
	//printf ("Moon yeclip = %.8f\n", yeclip);
	//printf ("Moon zeclip = %.8f\n", zeclip);

	double mlon =  rev(atan2d(yeclip, xeclip));
    //double latt  =  rev(atan2(zeclip, sqrt( xeclip*xeclip + yeclip*yeclip)));
    //r =  sqrt(xeclip*xeclip + yeclip*yeclip + zeclip*zeclip);
	//printf ("moon_long: Moon's longt is %f\n", mlon);

	//Compensate for Moon's perturbations
    //Sun's  mean longitude:        Ls     (already computed as Ls)
    //Moon's mean longitude:        Lm  =  N + w + Mm (for the Moon)
    //Sun's  mean anomaly:          Ms     (already computed as Ms)
    //Moon's mean anomaly:          Mm     (already computed in this function)
    //Moon's mean elongation:       D   =  Lm - Ls
    //Moon's argument of latitude:  F   =  Lm - N
	
	double Lm = rev(N + w + Mm);
	double D  = Lm - Ls;
	double F = Lm - N;
	//printf ("----------------\n");	
	//printf ("Sun's mean anomaly Ms = %.8f\n", Ms);
	//printf ("Moon mean anomaly Mm = %.8f\n", Mm);
	//printf ("Sun's mean longitude Ls = %.8f\n", Ls);
	//printf ("moon_long: Moon's mean longitude Lm is %f\n", Lm);
	//printf ("moon_long: Moon's argument of latitude F is %f\n", F);

	printf ("moon_long: Moon's longt before perturb fix is %f\n", mlon);
	mlon += 
		-1.274 * sind(Mm - 2*D)		//Evection
		+0.658 * sind(2*D)			//Variation
		-0.186 * sind(Ms)			//Yearly equation
		-0.059 * sind(2*Mm - 2*D)
		-0.057 * sind(Mm - 2*D + Ms)
		+0.053 * sind(Mm + 2*D)
		+0.046 * sind(2*D - Ms)
		+0.041 * sind(Mm - Ms)
		-0.035 * sind(D)			//Parallactic equation
		-0.031 * sind(Mm + Ms)
		-0.015 * sind(2*F - 2*D)
		+0.011 * sind(Mm - 4*D);
	mlon = rev(mlon);
	//printf ("Sun's longitude = %.8f\n", slon);
	printf ("moon_long: Moon's longt after perturb fix is %f\n", mlon);

	int n;

	//Calculate Tithi and Paksha
	if (mlon<slon)
		mlon += 360.0;
	printf ("Diff between Moons and Sun's longitudes = %.8f\n", mlon - slon);
	n = (int) (((mlon-slon)/12.0) + 0.0025);
	strcpy(pdata->tithi,__g_tithi[n]);
	printf ("Tithi index is n= %d -- %s\n", n, __g_tithi[n]);
	(n<=14)?strcpy(pdata->paksha,"Shukla"):strcpy(pdata->paksha,"Krishna");
}

void calculate_tithi_radians(int dd, int mm, int yy, double hr, double zhr, panchanga *pdata) {
	//https://www.stjarnhimlen.se/comp/tutorial.html#7
	double d = timescale (dd, mm, yy, hr, zhr);

	//------------------------------Sun
	//below are in radians
    double w = radians(rev(282.9404 + 4.70935e-5 * d));		//Sun's longitude of perihelion
    double Ms = radians(rev(356.0470 + 0.9856002585 * d));	//Sun's mean anomaly
	//double oblecl = 23.4393 - 3.563e-7 * d;		//obliquity of the ecliptic

    double e = 0.016709 - 1.151e-9 * d;		//eccentricity

	double E = Ms + e * sin(Ms) * (1 + e * cos(Ms));//eccentricity anomaly

	//Sun's mean longitude
	double Ls = w + Ms;
	//printf ("Sun's eccentricity anomaly E = %.8f\n", E);

	//Sun's rectangular coordinates
	double x = cos(E) - e;
    double y = sin(E) * sqrt(1 - e*e);

	//distance from Sun and true anomaly
	//double r = sqrt(x*x + y*y);	//in Earth radii
    double v = atan2( y, x );	//true anomaly
	double slon = v + w;
	printf ("Sun's longitude = %.8f\n", rev(degrees(slon)));

	//------------------------------Moon
	//all below are in radians
	double N = radians(rev(125.1228 - 0.0529538083  * d));   //Longt of ascending node
    const double i = 0.089804;			//Inclination in degrees is 5.1454
    w = radians(rev(318.0634 + 0.1643573223  * d));		//Arg. of perigee
    double Mm = radians(rev(115.3654 + 13.0649929509 * d));  //Mean eccentricity anomaly

    const double a = 60.2666; //Mean distance in Earth equatorial radii
    e = 0.054900;//Eccentricity

	//iterate for accurate eccentricity anomaly
	E = Mm + e * sin(Mm) * (1 + e * cos(Mm));
	double eps;
	do {
		eps	= (E - e * sin(E) - Mm) / (1 - e * cos(E));
		E = E - eps;
	} while (eps > 0.005 || eps < -0.005);
	//compute rectangular (x,y) coordinates in the plane of the lunar orbit
	x = a * (cos(E) - e);
    y = a * sqrt(1 - e*e) * sin(E);
	//printf ("Moon x = %.8f\n", x);
	//printf ("Moon y = %.8f\n", y);

	double r = sqrt(x*x + y*y); //distance Earth radii
    v = atan2(y, x); //true anomaly

	double xeclip = r * (cos(N) * cos(v+w) - sin(N) * sin(v+w) * cos(i));
    double yeclip = r * (sin(N) * cos(v+w) + cos(N) * sin(v+w) * cos(i));
    //double zeclip = r * sin(v+w) * sin(i);
	//printf ("Moon xeclip = %.8f\n", xeclip);
	//printf ("Moon yeclip = %.8f\n", yeclip);
	//printf ("Moon zeclip = %.8f\n", zeclip);

	double mlon =  atan2(yeclip, xeclip);
    //double latt  =  atan2(zeclip, sqrt( xeclip*xeclip + yeclip*yeclip));
    //r =  sqrt(xeclip*xeclip + yeclip*yeclip + zeclip*zeclip);

	//Compensate for Moon's perturbations
    //Sun's  mean longitude:        Ls     (already computed as Ls)
    //Moon's mean longitude:        Lm  =  N + w + Mm (for the Moon)
    //Sun's  mean anomaly:          Ms     (already computed as Ms)
    //Moon's mean anomaly:          Mm     (already computed in this function)
    //Moon's mean elongation:       D   =  Lm - Ls
    //Moon's argument of latitude:  F   =  Lm - N
	
	double Lm = N + w + Mm;
	double D  = Lm - Ls;
	double F = Lm - N;
	//printf ("----------------\n");	
	//printf ("Sun's mean anomaly Ms = %.8f\n", Ms);
	//printf ("Moon mean anomaly Mm = %.8f\n", Mm);
	//printf ("Sun's mean longitude Ls = %.8f\n", Ls);
	//printf ("moon_long: Moon's mean longitude Lm is %f\n", Lm);
	//printf ("moon_long: Moon's argument of latitude F is %f\n", F);

	printf ("moon_long: Moon's longt before perturb fix in radians is %f\n", mlon);
	printf ("moon_long: Moon's longt before perturb fix is %f\n", rev(degrees(mlon)));
	mlon += 
		-1.274 * sin(Mm - 2*D)		//Evection
		+0.658 * sin(2*D)			//Variation
		-0.186 * sin(Ms)			//Yearly equation
		-0.059 * sin(2*Mm - 2*D)
		-0.057 * sin(Mm - 2*D + Ms)
		+0.053 * sin(Mm + 2*D)
		+0.046 * sin(2*D - Ms)
		+0.041 * sin(Mm - Ms)
		-0.035 * sin(D)			//Parallactic equation
		-0.031 * sin(Mm + Ms)
		-0.015 * sin(2*F - 2*D)
		+0.011 * sin(Mm - 4*D);
	printf ("moon_long: Moon's longt after perturb fix in radians is %f\n", mlon);

	mlon = rev(degrees(mlon));
	//printf ("Sun's longitude = %.8f\n", slon);
	printf ("moon_long: Moon's longt after perturb fix is %f\n", mlon);

	int n;

	//Calculate Tithi and Paksha
	if (mlon<slon)
		mlon += 360.0;
	printf ("Diff between Moons and Sun's longitudes = %.8f\n", mlon - slon);
	n = (int) (((mlon-slon)/12.0) + 0.0025);
	strcpy(pdata->tithi,__g_tithi[n]);
	printf ("Tithi index is n= %d -- %s\n", n, __g_tithi[n]);
	(n<=14)?strcpy(pdata->paksha,"Shukla"):strcpy(pdata->paksha,"Krishna");
}

int main () {
	//https://www.satellite-calculations.com/Satellite/suncalc.htm
	//https://www.heavens-above.com/moon.aspx
	//https://www.stjarnhimlen.se/comp/tutorial.html
	
	panchanga pdata;	
	double srise, sset;
	double sunhrs = sun_esrl (&srise, &sset, 
					31, 10, 2019, 
					5.5, 12.98889, 77.62210, 0); //Bangalore

	//printf ("srise = %.15f\n", srise);
	//printf ("sset = %.15f\n", sset);
	printf ("srise = %02d:%02d\n", (int)floor(srise), (int)floor(fpart(srise)*60));
	printf ("sset = %02d:%02d\n", (int)floor(sset), (int)floor(fpart(sset)*60));
	//printf ("sunhrs = %.15f\n", sunhrs);
	//
	printf ("--------Alt--------\n");	
	sunhrs = sun_esrl (&srise, &sset, 
					31, 10, 2019, 
					5.5, 12.98889, 77.62210, 1); //Bangalore

	printf ("srise = %02d:%02d\n", (int)floor(srise), (int)floor(fpart(srise)*60));
	printf ("sset = %02d:%02d\n", (int)floor(sset), (int)floor(fpart(sset)*60));
	//printf ("sunhrs = %.15f\n", sunhrs);
	printf ("=======19, 4, 1990=========\n");	
	calculate_tithi (19, 4, 1990, 5.5, 5.5, &pdata);
	printf ("=======19, 4, 1990====radians=====\n");	
	calculate_tithi_radians (19, 4, 1990, 5.5, 5.5, &pdata);
	//printf ("======13, 10, 2019==========\n");	
	//calculate_tithi (13, 10, 2019, 0.6, 5.5, &pdata);
	//printf ("======13, 10, 2019==========\n");	
	//calculate_tithi (13, 10, 2019, 0.7, 5.5, &pdata);
	//printf ("======14, 10, 2019==========\n");	
	//calculate_tithi (14, 10, 2019, 2.6166667, 5.5, &pdata);
	//printf ("======15, 10, 2019==========\n");	
	//calculate_tithi (15, 10, 2019, 4.33333333, 5.5, &pdata);
	//printf ("======19, 10, 2019==========\n");	
	//calculate_tithi (19, 10, 2019, 7.7166667, 5.5, &pdata);
	//printf ("======20, 10, 2019==========\n");	
	//calculate_tithi (20, 10, 2019, 7.483333, 5.5, &pdata);
	//printf ("======20, 10, 2019==========\n");	
	//calculate_tithi (20, 10, 2019, 7.4666667, 5.5, &pdata);
	//printf ("======20, 10, 2019==========\n");	
	//calculate_tithi (20, 10, 2019, 7.45, 5.5, &pdata);
}

