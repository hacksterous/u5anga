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

double sun_esrl (double* sunriseLSTp, 
			double* sunsetLSTp, 
			int dd, 
			int mm, 
			int yyyy, 
			double zhr, 
			double latt, 
			double longt) {

	// Reverse engineered from the NOAA Excel:
	//https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html

	double julian_day;
	julian_day = JulianDay (dd, mm, yyyy, 12.0, zhr);
	//printf ("sun_esrl: julian_day is %.15f\n", julian_day);

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

double nround(double x, unsigned int n) {
    double fac = pow(10, n);
    return round(x*fac)/fac;
}

void calculate_tithi(int dd, int mm, int yy, double hr, double zhr, panchanga *pdata) {
	//https://www.stjarnhimlen.se/comp/tutorial.html#7
	//and also https://stjarnhimlen.se/comp/ppcomp.html
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
	double slon = (v + w);
	//printf ("Sun's longitude = %.8f\n", slon);

	//------------------------------Moon
	//all below are in degrees
	double N = rev(125.1228 - 0.0529538083  * d);   //Longt of ascending node
    const double i = 5.1454;						//Inclination
    w = rev(318.0634 + 0.1643573223  * d);			//Arg. of perigee
    double Mm = rev(115.3654 + 13.0649929509 * d);  //Mean eccentricity anomaly

    const double a = 60.2666; //Mean distance in Earth equatorial radii
    e = 0.054900;//Eccentricity

	//iterate for accurate eccentricity anomaly
	E = Mm + (180/PI) * e * sind(Mm) * (1 + e * cosd(Mm));
	double eps;
	do {
		eps	= (E - (180/PI) * e * sind(E) - Mm) / (1 - e * cosd(E));
		E = E - eps;
	} while (eps > 0.001 || eps < -0.001);
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

	double mlon =  (atan2d(yeclip, xeclip));
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

	//printf ("moon_long: Moon's longt before perturb fix is %f\n", mlon);
	mlon += 
		-1.274 * sind(Mm - 2*D)		//Evection
		+0.658 * sind(2*D)			//Variation
		-0.186 * sind(Ms)			//Yearly equation
		-0.059 * sind(2*Mm - 2*D)
		-0.057 * sind(Mm - 2*D + Ms)
		+0.053 * sind(Mm + 2*D)
		+0.046 * sind(2*D - Ms)
		+0.041 * sind(Mm - Ms)
		-0.035 * sind(D)			//Parallactic equation -- should be -0.0375
		-0.0025 * sind(D)			//Parallactic equation -- new term
									//anirb's term actual val from https://en.wikipedia.org/wiki/Lunar_theory
									//is 125" = 0.0375 degrees
		-0.031 * sind(Mm + Ms)
		-0.0144 * sind(2*F - 2*D)   //reduction to the ecliptic was 0.015
		+0.011 * sind(Mm - 4*D);
	//printf ("Sun's longitude = %.8f\n", slon);
	//printf ("moon_long: Moon's longt after perturb fix is %f\n", mlon);

	int n;

	//Calculate Tithi and Paksha
	//printf ("Diff between Moons and Sun's longitudes = %.8f\n", mlon - slon);
	//printf ("Index of diff between Moons and Sun's longitudes = %.8f\n", (mlon - slon)/12.0);
	if (mlon<slon) {
		mlon += 360.0;
		n = (int) (nround((rev(mlon-slon+0.07)/12.0), 3));
	} else
		n = (int) (nround((rev(mlon-slon)/12.0), 3));
	strcpy(pdata->tithi,__g_tithi[n]);
	//printf ("Tithi index is n= %d -- %s\n", n, __g_tithi[n]);
	printf ("%s\n", __g_tithi[n]);
	//printf ("%f\n", nround(rev(mlon-slon+0.06)/12.0, 3));
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
    double v = atan2( y, x );		//true anomaly
	double slon = v + w;
	//printf ("Sun's longitude = %.8f\n", rev(degrees(slon)));

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
	} while (eps > 0.001 || eps < -0.001);
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

	//printf ("moon_long: Moon's longt before perturb fix in radians is %f\n", mlon);
	//printf ("moon_long: Moon's longt before perturb fix is %f\n", rev(degrees(mlon)));
	mlon += //in radians
		-0.022235495 * sin(Mm - 2*D)	//Evection
		+0.011484266 * sin(2*D)			//Variation
		-0.003246312 * sin(Ms)			//Yearly equation
		-0.001029744 * sin(2*Mm - 2*D)
		-0.000994838 * sin(Mm - 2*D + Ms)
		+0.000925025 * sin(Mm + 2*D)
		+0.000802851 * sin(2*D - Ms)
		+0.000715585 * sin(Mm - Ms)
		-0.000610865 * sin(D)			//Parallactic equation
		-4.3633231292e-5 * sin(D)		//Parallactic equation -- new term
										//anirb's term actual val from https://en.wikipedia.org/wiki/Lunar_theory
										//is 125" = 0.0375 degrees
		-0.000541052 * sin(Mm + Ms)
		-0.0002513274 * sin(2*F - 2*D)	//reduction to the ecliptic was 0.000261799    
		+0.000191986 * sin(Mm - 4*D);
	//printf ("moon_long: Moon's longt after perturb fix in radians is %f\n", mlon);

	slon = (degrees(slon));
	mlon = (degrees(mlon));
	//printf ("Sun's longitude = %.8f\n", slon);
	//printf ("moon_long: Moon's longt after perturb fix is %f\n", mlon);

	int n;

	//Calculate Tithi and Paksha
	if (mlon<slon) {
		mlon += 360.0;
		n = (int) (nround((rev(mlon-slon+0.07)/12.0), 3));
	} else
		n = (int) (nround((rev(mlon-slon)/12.0), 3));
	//printf ("Diff between Moons and Sun's longitudes = %.8f\n", mlon - slon);
	//printf ("Index of diff between Moons and Sun's longitudes = %.8f\n", (mlon - slon)/12.0);
	strcpy(pdata->tithi,__g_tithi[n]);
	//printf ("Tithi index is n= %d -- %s\n", n, __g_tithi[n]);
	printf ("%s\n", __g_tithi[n]);
	(n<=14)?strcpy(pdata->paksha,"Shukla"):strcpy(pdata->paksha,"Krishna");
}

int main () {
	//https://www.satellite-calculations.com/Satellite/suncalc.htm
	//https://www.heavens-above.com/moon.aspx
	//https://www.stjarnhimlen.se/comp/tutorial.html
	
	panchanga pdata;	
	double srise, sset, sunhrs;
	//printf ("-------------------\n");	
	sunhrs = sun_esrl (&srise, &sset, 
					31, 10, 2019, 
					5.5, 12.98889, 77.62210); //Bangalore

	//printf ("srise = %02d:%02d\n", (int)floor(srise), (int)floor(fpart(srise)*60));
	//printf ("sset = %02d:%02d\n", (int)floor(sset), (int)floor(fpart(sset)*60));
	//printf ("-------------------\n");	
	//printf ("sunhrs = %.15f\n", sunhrs);
	//printf ("=======19, 4, 1990=========\n");	
	//calculate_tithi (19, 4, 1990, 5.5, 5.5, &pdata);
	//printf ("=======19, 4, 1990====radians=====\n");	
	//calculate_tithi_radians (19, 4, 1990, 5.5, 5.5, &pdata);

calculate_tithi_radians(   1 , 1, 2004, 6.35, 5.5, &pdata);
calculate_tithi_radians(   2 , 1, 2004, 8.866666666666667, 5.5, &pdata);
calculate_tithi_radians(   3 , 1, 2004, 11.616666666666667, 5.5, &pdata);
calculate_tithi_radians(   4 , 1, 2004, 14.383333333333333, 5.5, &pdata);
calculate_tithi_radians(   5 , 1, 2004, 16.966666666666665, 5.5, &pdata);
calculate_tithi_radians(   6 , 1, 2004, 19.25, 5.5, &pdata);
calculate_tithi_radians(   7 , 1, 2004, 21.166666666666668, 5.5, &pdata);
calculate_tithi_radians(   8 , 1, 2004, 22.65, 5.5, &pdata);
calculate_tithi_radians(   9 , 1, 2004, 23.716666666666665, 5.5, &pdata);
calculate_tithi_radians(   11, 1, 2004, 0.38333333333333336, 5.5, &pdata);
calculate_tithi_radians(   12, 1, 2004, 0.6333333333333333, 5.5, &pdata);
calculate_tithi_radians(   13, 1, 2004, 0.48333333333333334, 5.5, &pdata);
calculate_tithi_radians(   13, 1, 2004, 23.916666666666668, 5.5, &pdata);
calculate_tithi_radians(   14, 1, 2004, 22.916666666666668, 5.5, &pdata);
calculate_tithi_radians(   15, 1, 2004, 21.466666666666665, 5.5, &pdata);
calculate_tithi_radians(   16, 1, 2004, 19.583333333333332, 5.5, &pdata);
calculate_tithi_radians(   17, 1, 2004, 17.283333333333335, 5.5, &pdata);
calculate_tithi_radians(   18, 1, 2004, 14.633333333333333, 5.5, &pdata);
calculate_tithi_radians(   19, 1, 2004, 11.7, 5.5, &pdata);
calculate_tithi_radians(   20, 1, 2004, 8.616666666666667, 5.5, &pdata);
calculate_tithi_radians(   21, 1, 2004, 5.516666666666667, 5.5, &pdata);
calculate_tithi_radians(   22, 1, 2004, 2.5666666666666664, 5.5, &pdata);
calculate_tithi_radians(   22, 1, 2004, 23.95, 5.5, &pdata);
calculate_tithi_radians(   23, 1, 2004, 21.8, 5.5, &pdata);
calculate_tithi_radians(   24, 1, 2004, 20.316666666666666, 5.5, &pdata);
calculate_tithi_radians(   25, 1, 2004, 19.6, 5.5, &pdata);
calculate_tithi_radians(   26, 1, 2004, 19.733333333333334, 5.5, &pdata);
calculate_tithi_radians(   27, 1, 2004, 20.7, 5.5, &pdata);
calculate_tithi_radians(   28, 1, 2004, 22.45, 5.5, &pdata);
calculate_tithi_radians(   30, 1, 2004, 0.7666666666666667, 5.5, &pdata);
calculate_tithi_radians(   31, 1, 2004, 3.4333333333333336, 5.5, &pdata);
calculate_tithi_radians(   1 , 2, 2004, 6.166666666666667, 5.5, &pdata);
calculate_tithi_radians(   2 , 2, 2004, 8.716666666666667, 5.5, &pdata);
calculate_tithi_radians(   3 , 2, 2004, 10.883333333333333, 5.5, &pdata);
calculate_tithi_radians(   4 , 2, 2004, 12.566666666666666, 5.5, &pdata);
calculate_tithi_radians(   5 , 2, 2004, 13.683333333333334, 5.5, &pdata);
calculate_tithi_radians(   6 , 2, 2004, 14.266666666666667, 5.5, &pdata);
calculate_tithi_radians(   7 , 2, 2004, 14.366666666666667, 5.5, &pdata);
calculate_tithi_radians(   8 , 2, 2004, 14.05, 5.5, &pdata);
calculate_tithi_radians(   9 , 2, 2004, 13.383333333333333, 5.5, &pdata);
calculate_tithi_radians(   10, 2, 2004, 12.416666666666666, 5.5, &pdata);
calculate_tithi_radians(   11, 2, 2004, 11.2, 5.5, &pdata);
calculate_tithi_radians(   12, 2, 2004, 9.75, 5.5, &pdata);
calculate_tithi_radians(   13, 2, 2004, 8.066666666666666, 5.5, &pdata);
calculate_tithi_radians(   14, 2, 2004, 6.183333333333334, 5.5, &pdata);
calculate_tithi_radians(   15, 2, 2004, 4.083333333333333, 5.5, &pdata);
calculate_tithi_radians(   16, 2, 2004, 1.8333333333333335, 5.5, &pdata);
calculate_tithi_radians(   16, 2, 2004, 23.45, 5.5, &pdata);
calculate_tithi_radians(   17, 2, 2004, 21.05, 5.5, &pdata);
calculate_tithi_radians(   18, 2, 2004, 18.716666666666665, 5.5, &pdata);
calculate_tithi_radians(   19, 2, 2004, 16.583333333333332, 5.5, &pdata);
calculate_tithi_radians(   20, 2, 2004, 14.783333333333333, 5.5, &pdata);
calculate_tithi_radians(   21, 2, 2004, 13.466666666666667, 5.5, &pdata);
calculate_tithi_radians(   22, 2, 2004, 12.766666666666667, 5.5, &pdata);
calculate_tithi_radians(   23, 2, 2004, 12.766666666666667, 5.5, &pdata);
calculate_tithi_radians(   24, 2, 2004, 13.516666666666667, 5.5, &pdata);
calculate_tithi_radians(   25, 2, 2004, 14.983333333333333, 5.5, &pdata);
calculate_tithi_radians(   26, 2, 2004, 17.066666666666666, 5.5, &pdata);
calculate_tithi_radians(   27, 2, 2004, 19.566666666666666, 5.5, &pdata);
calculate_tithi_radians(   28, 2, 2004, 22.216666666666665, 5.5, &pdata);
calculate_tithi_radians(   1 , 3, 2004, 0.7666666666666667, 5.5, &pdata);
calculate_tithi_radians(   2 , 3, 2004, 2.95, 5.5, &pdata);
calculate_tithi_radians(   3 , 3, 2004, 4.583333333333333, 5.5, &pdata);
calculate_tithi_radians(   4 , 3, 2004, 5.566666666666666, 5.5, &pdata);
calculate_tithi_radians(   5 , 3, 2004, 5.9, 5.5, &pdata);
calculate_tithi_radians(   6 , 3, 2004, 5.583333333333333, 5.5, &pdata);
calculate_tithi_radians(   7 , 3, 2004, 4.733333333333333, 5.5, &pdata);
calculate_tithi_radians(   8 , 3, 2004, 3.4333333333333336, 5.5, &pdata);
calculate_tithi_radians(   9 , 3, 2004, 1.7833333333333332, 5.5, &pdata);
calculate_tithi_radians(   9 , 3, 2004, 23.9, 5.5, &pdata);
calculate_tithi_radians(   10, 3, 2004, 21.883333333333333, 5.5, &pdata);
calculate_tithi_radians(   11, 3, 2004, 19.783333333333335, 5.5, &pdata);
calculate_tithi_radians(   12, 3, 2004, 17.65, 5.5, &pdata);
calculate_tithi_radians(   13, 3, 2004, 15.55, 5.5, &pdata);
calculate_tithi_radians(   14, 3, 2004, 13.483333333333333, 5.5, &pdata);
calculate_tithi_radians(   15, 3, 2004, 11.5, 5.5, &pdata);
calculate_tithi_radians(   16, 3, 2004, 9.633333333333333, 5.5, &pdata);
calculate_tithi_radians(   17, 3, 2004, 7.916666666666667, 5.5, &pdata);
calculate_tithi_radians(   18, 3, 2004, 6.45, 5.5, &pdata);
calculate_tithi_radians(   19, 3, 2004, 5.266666666666667, 5.5, &pdata);
calculate_tithi_radians(   20, 3, 2004, 4.483333333333333, 5.5, &pdata);
calculate_tithi_radians(   21, 3, 2004, 4.183333333333334, 5.5, &pdata);
calculate_tithi_radians(   22, 3, 2004, 4.416666666666667, 5.5, &pdata);
calculate_tithi_radians(   23, 3, 2004, 5.25, 5.5, &pdata);
calculate_tithi_radians(   24, 3, 2004, 6.666666666666667, 5.5, &pdata);
calculate_tithi_radians(   25, 3, 2004, 8.616666666666667, 5.5, &pdata);
calculate_tithi_radians(   26, 3, 2004, 10.983333333333333, 5.5, &pdata);
calculate_tithi_radians(   27, 3, 2004, 13.55, 5.5, &pdata);
calculate_tithi_radians(   28, 3, 2004, 16.1, 5.5, &pdata);
calculate_tithi_radians(   29, 3, 2004, 18.383333333333333, 5.5, &pdata);
calculate_tithi_radians(   30, 3, 2004, 20.166666666666668, 5.5, &pdata);
calculate_tithi_radians(   31, 3, 2004, 21.316666666666666, 5.5, &pdata);
calculate_tithi_radians(   1 , 4, 2004, 21.716666666666665, 5.5, &pdata);
calculate_tithi_radians(   2 , 4, 2004, 21.366666666666667, 5.5, &pdata);
calculate_tithi_radians(   3 , 4, 2004, 20.333333333333332, 5.5, &pdata);
calculate_tithi_radians(   4 , 4, 2004, 18.683333333333334, 5.5, &pdata);
calculate_tithi_radians(   5 , 4, 2004, 16.533333333333335, 5.5, &pdata);
calculate_tithi_radians(   6 , 4, 2004, 14.033333333333333, 5.5, &pdata);
calculate_tithi_radians(   7 , 4, 2004, 11.3, 5.5, &pdata);
calculate_tithi_radians(   8 , 4, 2004, 8.483333333333333, 5.5, &pdata);
calculate_tithi_radians(   9 , 4, 2004, 5.666666666666667, 5.5, &pdata);
calculate_tithi_radians(   10, 4, 2004, 2.966666666666667, 5.5, &pdata);
calculate_tithi_radians(   11, 4, 2004, 0.4666666666666667, 5.5, &pdata);
calculate_tithi_radians(   11, 4, 2004, 22.25, 5.5, &pdata);
calculate_tithi_radians(   12, 4, 2004, 20.35, 5.5, &pdata);
calculate_tithi_radians(   13, 4, 2004, 18.833333333333332, 5.5, &pdata);
calculate_tithi_radians(   14, 4, 2004, 17.716666666666665, 5.5, &pdata);
calculate_tithi_radians(   15, 4, 2004, 17.033333333333335, 5.5, &pdata);
calculate_tithi_radians(   16, 4, 2004, 16.783333333333335, 5.5, &pdata);
calculate_tithi_radians(   17, 4, 2004, 17.0, 5.5, &pdata);
calculate_tithi_radians(   18, 4, 2004, 17.683333333333334, 5.5, &pdata);
calculate_tithi_radians(   19, 4, 2004, 18.85, 5.5, &pdata);
calculate_tithi_radians(   20, 4, 2004, 20.45, 5.5, &pdata);
calculate_tithi_radians(   21, 4, 2004, 22.45, 5.5, &pdata);
calculate_tithi_radians(   23, 4, 2004, 0.7833333333333333, 5.5, &pdata);
calculate_tithi_radians(   24, 4, 2004, 3.3, 5.5, &pdata);
calculate_tithi_radians(   25, 4, 2004, 5.85, 5.5, &pdata);
calculate_tithi_radians(   26, 4, 2004, 8.233333333333333, 5.5, &pdata);
calculate_tithi_radians(   27, 4, 2004, 10.233333333333333, 5.5, &pdata);
calculate_tithi_radians(   28, 4, 2004, 11.683333333333334, 5.5, &pdata);
calculate_tithi_radians(   29, 4, 2004, 12.433333333333334, 5.5, &pdata);
calculate_tithi_radians(   30, 4, 2004, 12.433333333333334, 5.5, &pdata);
calculate_tithi_radians(   1 , 5, 2004, 11.633333333333333, 5.5, &pdata);
calculate_tithi_radians(   2 , 5, 2004, 10.1, 5.5, &pdata);
calculate_tithi_radians(   3 , 5, 2004, 7.9, 5.5, &pdata);
calculate_tithi_radians(   4 , 5, 2004, 5.183333333333334, 5.5, &pdata);
calculate_tithi_radians(   5 , 5, 2004, 2.05, 5.5, &pdata);
calculate_tithi_radians(   5 , 5, 2004, 22.666666666666668, 5.5, &pdata);
calculate_tithi_radians(   6 , 5, 2004, 19.2, 5.5, &pdata);
calculate_tithi_radians(   7 , 5, 2004, 15.766666666666667, 5.5, &pdata);
calculate_tithi_radians(   8 , 5, 2004, 12.55, 5.5, &pdata);
calculate_tithi_radians(   9 , 5, 2004, 9.65, 5.5, &pdata);
calculate_tithi_radians(   10, 5, 2004, 7.2, 5.5, &pdata);
calculate_tithi_radians(   11, 5, 2004, 5.3, 5.5, &pdata);
calculate_tithi_radians(   12, 5, 2004, 3.9833333333333334, 5.5, &pdata);
calculate_tithi_radians(   13, 5, 2004, 3.283333333333333, 5.5, &pdata);
calculate_tithi_radians(   14, 5, 2004, 3.216666666666667, 5.5, &pdata);
calculate_tithi_radians(   15, 5, 2004, 3.7333333333333334, 5.5, &pdata);
calculate_tithi_radians(   16, 5, 2004, 4.783333333333333, 5.5, &pdata);
calculate_tithi_radians(   17, 5, 2004, 6.283333333333333, 5.5, &pdata);
calculate_tithi_radians(   18, 5, 2004, 8.166666666666666, 5.5, &pdata);

printf ("-------------------\n");	

}

