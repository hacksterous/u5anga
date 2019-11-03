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

static char __g_tithi[][15]={"Pratipad","Dvitiya","Tritiya","Chaturthi","Panchami",
		  "Shashthi","Saptami","Ashtami","Navami","Dashami","Ekadashi",
		  "Dvadashi","Trayodashi","Chaturdashi","Purnima","Pratipad",
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
				double hr, 
				double zhr) {
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
	return (double) (d + (hr - zhr)/24.0);//add 2451543.5 to get Julian Date;
}

double radians (double degrees) {
	return (3.14159265358979323846 * degrees / 180.0);
}

double degrees (double radians) {
	return (180.0 * radians /3.14159265358979323846);
}

double fpart (double x) {
	double y = x - floor(x);
	if (y < 0)
		y += 1;
	return y;
}

double JulianDay (int date, int month, int year, double hr, double zhr) {
    if (month<=2) {month=month+12; year=year-1;}

	//Use formula given in 
	//A Physically-Based Night Sky ModelA Physically-Based Night Sky Model
	//by Julie Dorsey
	
	//JD = 1720996.5 − floor(year/100) + floor(year/400) + 
	//		floor(365.25 * year) + floor (30.6001 * (month + 1)) + 
	//		D + (h + (m + s/60)/60)/24
	//		### h:m:s is local time with zone correction
	//NOTE: deltaT of 65 seconds (=0.018055555555556hr) has not been added to hr

    return floor(year/400) - floor (year/100) + floor(365.25*year) + 
			floor(30.6001*(month+1)) + 1720996.5 + date + (hr - zhr)/24.0;
}

double sun_esrl (double* sunriseLSTp, 
			double* sunsetLSTp, 
			double julian_day, 
			double zhr, 
			double latt, 
			double longt) {

	//Reverse engineered from the NOAA Excel:
	//https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html

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
	//printf ("sun_true_long_deg is %.15f\n", rev(sun_true_long_deg));

	double sun_true_anom_deg = sun_geom_mean_anom_deg + sun_eqn_of_centre; //N
	//printf ("sun_true_anom_deg is %.15f\n", sun_true_anom_deg);

	double sun_rad_vector_au = (1.000001018*(1-earth_orbit_eccentricity*earth_orbit_eccentricity))/
							(1+earth_orbit_eccentricity*cos(radians(sun_true_anom_deg))); //O
	//printf ("sun_rad_vector_au is %.15f\n", sun_rad_vector_au);

	double sun_apparent_long_deg = sun_true_long_deg-0.00569-0.00478*
							sin(radians(125.04-1934.136*julian_century)); //P
	//printf ("sun_apparent_long_deg is %.15f\n", rev(sun_apparent_long_deg));

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

	double hour_angle_sunrise_deg = degrees(acos(cos(radians(90.833))/(cos(radians(latt))*
							cos(radians(sun_declination_deg)))-tan(radians(latt))*
							tan(radians(sun_declination_deg)))); //W
	//printf ("hour_angle_sunrise_deg is %.15f\n", hour_angle_sunrise_deg);

	double solar_noon_LST = (720.0-4*longt-equation_of_time_min+zhr*60)/1440.0; //X
	//printf ("solar_noon_LST is %.15f\n", solar_noon_LST);

	*sunriseLSTp = (solar_noon_LST - (hour_angle_sunrise_deg*4)/1440.0) * 24; //Y

	*sunsetLSTp = (solar_noon_LST + (hour_angle_sunrise_deg*4)/1440.0)*24; //Z

	double sunlight_duration_hrs = hour_angle_sunrise_deg * 8 / 60; //AA

	return (sunlight_duration_hrs);
}

double nround(double x, unsigned int fac) {
    return round(x*fac)/fac;
}

void calculate_tithi_radians(int dd, int mm, int yy, double hr, double zhr, panchanga *pdata) {
	//https://www.stjarnhimlen.se/comp/tutorial.html#7
	double d = timescale (dd, mm, yy, hr, zhr);
	printf ("calculate_tithi_radians d= %.8f\n", d);

	//------------------------------Sun
	//below are in radians
    double w = radians(rev(282.9404 + 4.70935e-5 * d));		//Sun's longitude of perihelion
    double Ms = radians(rev(356.0470 + 0.9856002585 * d));	//Sun's mean anomaly
	//double oblecl = 23.4393 - 3.563e-7 * d;		//obliquity of the ecliptic

    double e = 0.016709 - 1.151e-9 * d;		//eccentricity

	double E = Ms + e * sin(Ms) * (1 + e * cos(Ms));//eccentricity anomaly

	//Sun's mean longitude
	double Ls = w + Ms;

	//Sun's rectangular coordinates
	double x = cos(E) - e;
    double y = sin(E) * sqrt(1 - e*e);

	//distance from Sun and true anomaly
	//double r = sqrt(x*x + y*y);	//in Earth radii
    double v = atan2( y, x );		//true anomaly
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
	} while (eps > 0.001 || eps < -0.001);
	//compute rectangular (x,y) coordinates in the plane of the lunar orbit
	x = a * (cos(E) - e);
    y = a * sqrt(1 - e*e) * sin(E);

	double r = sqrt(x*x + y*y); //distance Earth radii
    v = atan2(y, x); //true anomaly

	double xeclip = r * (cos(N) * cos(v+w) - sin(N) * sin(v+w) * cos(i));
    double yeclip = r * (sin(N) * cos(v+w) + cos(N) * sin(v+w) * cos(i));
    //double zeclip = r * sin(v+w) * sin(i);

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
	printf ("moon_long: Moon's longt before perturb fix is %f\n", rev(degrees(mlon)));
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
		-0.0002513274 * sin(2*F - 2*D)	//reduction to the ecliptic -- was 0.000261799    
		+0.000191986 * sin(Mm - 4*D);
	//printf ("moon_long: Moon's longt after perturb fix in radians is %f\n", mlon);

	slon = rev(degrees(slon));
	mlon = rev(degrees(mlon));

	printf ("Sun's ecl. longitude = %.8f\n", slon);
	printf ("Moon's ecl. longitude = %.8f\n", mlon);

	printf ("Sun's longitude = %.8f\n", slon);
	printf ("moon_long: Moon's longt after perturb fix is %f\n", mlon);

	int n;

	//Calculate Tithi and Paksha
	n = (int) (nround((rev(mlon-slon)/12.0), 1000));
	printf ("Diff between Moons and Sun's longitudes = %.8f\n", mlon - slon);
	printf ("Index of diff between Moons and Sun's longitudes = %.8f and index = %d\n", rev(mlon - slon)/12.0, n);
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
	double srise, sset, sunlightduration;

	double julian_day;
	julian_day = JulianDay (15, 10, 2019, 12.0, 5.5);
	//printf ("sun_esrl: julianday julian_day is %.15f\n", julian_day);

	julian_day =  timescale (15, 10, 2019, 12.0, 5.5) + 2451543.5;
	//printf ("sun_esrl: timescale julian_day is %.15f\n", julian_day);

	sunlightduration = sun_esrl (&srise, &sset, 
					julian_day,
					5.5, 12.93340, 77.59630); //Bangalore
	//printf ("srise -- %f\n", srise);
	//printf ("srise = %02d:%02d:%02d\n", (int)srise, (int)(fpart(srise)*60), (int)nround(fpart(fpart(srise)*60)*60, 1000));
	//printf ("sset = %02d:%02d:%02d\n", (int)sset, (int)(fpart(sset)*60), (int)nround(fpart(fpart(sset)*60)*60, 1000));
	//printf ("-------------------\n");	
	//printf ("sunlightduration = %02d:%02d:%02d\n", (int)(sunlightduration/60), (int)(fpart(sunlightduration/60)*60), (int)nround(fpart(fpart(sunlightduration/60)*60)*60, 1000));
	//printf ("-------------------\n");	
	//calculate_tithi_radians (14, 10, 2019, 2.566667, 5.5, &pdata);	//02:34
	//calculate_tithi_radians (14, 10, 2019, 2.58333, 5.5, &pdata);	//02:35
	//calculate_tithi_radians (14, 10, 2019, 2.6, 5.5, &pdata);		//02:36
	calculate_tithi_radians (21, 10, 2019, 6.733333, 5.5, &pdata);
	printf ("-------------------\n");	

//calculate_tithi_radians(   1 , 1, 2004, 6.35, 5.5, &pdata);
//calculate_tithi_radians(   2 , 1, 2004, 8.866666666666667, 5.5, &pdata);
//calculate_tithi_radians(   3 , 1, 2004, 11.616666666666667, 5.5, &pdata);
//calculate_tithi_radians(   4 , 1, 2004, 14.383333333333333, 5.5, &pdata);
//calculate_tithi_radians(   5 , 1, 2004, 16.966666666666665, 5.5, &pdata);
//calculate_tithi_radians(   6 , 1, 2004, 19.25, 5.5, &pdata);
//calculate_tithi_radians(   7 , 1, 2004, 21.166666666666668, 5.5, &pdata);
//calculate_tithi_radians(   8 , 1, 2004, 22.65, 5.5, &pdata);
//calculate_tithi_radians(   9 , 1, 2004, 23.716666666666665, 5.5, &pdata);
//calculate_tithi_radians(   11, 1, 2004, 0.38333333333333336, 5.5, &pdata);
//calculate_tithi_radians(   12, 1, 2004, 0.6333333333333333, 5.5, &pdata);
//calculate_tithi_radians(   13, 1, 2004, 0.48333333333333334, 5.5, &pdata);
//calculate_tithi_radians(   13, 1, 2004, 23.916666666666668, 5.5, &pdata);
//calculate_tithi_radians(   14, 1, 2004, 22.916666666666668, 5.5, &pdata);
//calculate_tithi_radians(   15, 1, 2004, 21.466666666666665, 5.5, &pdata);
//calculate_tithi_radians(   16, 1, 2004, 19.583333333333332, 5.5, &pdata);
//calculate_tithi_radians(   17, 1, 2004, 17.283333333333335, 5.5, &pdata);
//calculate_tithi_radians(   18, 1, 2004, 14.633333333333333, 5.5, &pdata);
//calculate_tithi_radians(   19, 1, 2004, 11.7, 5.5, &pdata);
//calculate_tithi_radians(   20, 1, 2004, 8.616666666666667, 5.5, &pdata);
//calculate_tithi_radians(   21, 1, 2004, 5.516666666666667, 5.5, &pdata);
//calculate_tithi_radians(   22, 1, 2004, 2.5666666666666664, 5.5, &pdata);
//calculate_tithi_radians(   22, 1, 2004, 23.95, 5.5, &pdata);
//calculate_tithi_radians(   23, 1, 2004, 21.8, 5.5, &pdata);
//calculate_tithi_radians(   24, 1, 2004, 20.316666666666666, 5.5, &pdata);
//calculate_tithi_radians(   25, 1, 2004, 19.6, 5.5, &pdata);
//calculate_tithi_radians(   26, 1, 2004, 19.733333333333334, 5.5, &pdata);
//calculate_tithi_radians(   27, 1, 2004, 20.7, 5.5, &pdata);
//calculate_tithi_radians(   28, 1, 2004, 22.45, 5.5, &pdata);
//calculate_tithi_radians(   30, 1, 2004, 0.7666666666666667, 5.5, &pdata);
//calculate_tithi_radians(   31, 1, 2004, 3.4333333333333336, 5.5, &pdata);
//calculate_tithi_radians(   1 , 2, 2004, 6.166666666666667, 5.5, &pdata);
//calculate_tithi_radians(   2 , 2, 2004, 8.716666666666667, 5.5, &pdata);
//calculate_tithi_radians(   3 , 2, 2004, 10.883333333333333, 5.5, &pdata);
//calculate_tithi_radians(   4 , 2, 2004, 12.566666666666666, 5.5, &pdata);
//calculate_tithi_radians(   5 , 2, 2004, 13.683333333333334, 5.5, &pdata);
//calculate_tithi_radians(   6 , 2, 2004, 14.266666666666667, 5.5, &pdata);
//calculate_tithi_radians(   7 , 2, 2004, 14.366666666666667, 5.5, &pdata);
//calculate_tithi_radians(   8 , 2, 2004, 14.05, 5.5, &pdata);
//calculate_tithi_radians(   9 , 2, 2004, 13.383333333333333, 5.5, &pdata);
//calculate_tithi_radians(   10, 2, 2004, 12.416666666666666, 5.5, &pdata);
//calculate_tithi_radians(   11, 2, 2004, 11.2, 5.5, &pdata);
//calculate_tithi_radians(   12, 2, 2004, 9.75, 5.5, &pdata);
//calculate_tithi_radians(   13, 2, 2004, 8.066666666666666, 5.5, &pdata);
//calculate_tithi_radians(   14, 2, 2004, 6.183333333333334, 5.5, &pdata);
//calculate_tithi_radians(   15, 2, 2004, 4.083333333333333, 5.5, &pdata);
//calculate_tithi_radians(   16, 2, 2004, 1.8333333333333335, 5.5, &pdata);
//calculate_tithi_radians(   16, 2, 2004, 23.45, 5.5, &pdata);
//calculate_tithi_radians(   17, 2, 2004, 21.05, 5.5, &pdata);
//calculate_tithi_radians(   18, 2, 2004, 18.716666666666665, 5.5, &pdata);
//calculate_tithi_radians(   19, 2, 2004, 16.583333333333332, 5.5, &pdata);
//calculate_tithi_radians(   20, 2, 2004, 14.783333333333333, 5.5, &pdata);
//calculate_tithi_radians(   21, 2, 2004, 13.466666666666667, 5.5, &pdata);
//calculate_tithi_radians(   22, 2, 2004, 12.766666666666667, 5.5, &pdata);
//calculate_tithi_radians(   23, 2, 2004, 12.766666666666667, 5.5, &pdata);
//calculate_tithi_radians(   24, 2, 2004, 13.516666666666667, 5.5, &pdata);
//calculate_tithi_radians(   25, 2, 2004, 14.983333333333333, 5.5, &pdata);
//calculate_tithi_radians(   26, 2, 2004, 17.066666666666666, 5.5, &pdata);
//calculate_tithi_radians(   27, 2, 2004, 19.566666666666666, 5.5, &pdata);
//calculate_tithi_radians(   28, 2, 2004, 22.216666666666665, 5.5, &pdata);
//calculate_tithi_radians(   1 , 3, 2004, 0.7666666666666667, 5.5, &pdata);
//calculate_tithi_radians(   2 , 3, 2004, 2.95, 5.5, &pdata);
//calculate_tithi_radians(   3 , 3, 2004, 4.583333333333333, 5.5, &pdata);
//calculate_tithi_radians(   4 , 3, 2004, 5.566666666666666, 5.5, &pdata);
//calculate_tithi_radians(   5 , 3, 2004, 5.9, 5.5, &pdata);
//calculate_tithi_radians(   6 , 3, 2004, 5.583333333333333, 5.5, &pdata);
//calculate_tithi_radians(   7 , 3, 2004, 4.733333333333333, 5.5, &pdata);
//calculate_tithi_radians(   8 , 3, 2004, 3.4333333333333336, 5.5, &pdata);
//calculate_tithi_radians(   9 , 3, 2004, 1.7833333333333332, 5.5, &pdata);
//calculate_tithi_radians(   9 , 3, 2004, 23.9, 5.5, &pdata);
//calculate_tithi_radians(   10, 3, 2004, 21.883333333333333, 5.5, &pdata);
//calculate_tithi_radians(   11, 3, 2004, 19.783333333333335, 5.5, &pdata);
//calculate_tithi_radians(   12, 3, 2004, 17.65, 5.5, &pdata);
//calculate_tithi_radians(   13, 3, 2004, 15.55, 5.5, &pdata);
//calculate_tithi_radians(   14, 3, 2004, 13.483333333333333, 5.5, &pdata);
//calculate_tithi_radians(   15, 3, 2004, 11.5, 5.5, &pdata);
//calculate_tithi_radians(   16, 3, 2004, 9.633333333333333, 5.5, &pdata);
//calculate_tithi_radians(   17, 3, 2004, 7.916666666666667, 5.5, &pdata);
//calculate_tithi_radians(   18, 3, 2004, 6.45, 5.5, &pdata);
//calculate_tithi_radians(   19, 3, 2004, 5.266666666666667, 5.5, &pdata);
//calculate_tithi_radians(   20, 3, 2004, 4.483333333333333, 5.5, &pdata);
//calculate_tithi_radians(   21, 3, 2004, 4.183333333333334, 5.5, &pdata);
//calculate_tithi_radians(   22, 3, 2004, 4.416666666666667, 5.5, &pdata);
//calculate_tithi_radians(   23, 3, 2004, 5.25, 5.5, &pdata);
//calculate_tithi_radians(   24, 3, 2004, 6.666666666666667, 5.5, &pdata);
//calculate_tithi_radians(   25, 3, 2004, 8.616666666666667, 5.5, &pdata);
//calculate_tithi_radians(   26, 3, 2004, 10.983333333333333, 5.5, &pdata);
//calculate_tithi_radians(   27, 3, 2004, 13.55, 5.5, &pdata);
//calculate_tithi_radians(   28, 3, 2004, 16.1, 5.5, &pdata);
//calculate_tithi_radians(   29, 3, 2004, 18.383333333333333, 5.5, &pdata);
//calculate_tithi_radians(   30, 3, 2004, 20.166666666666668, 5.5, &pdata);
//calculate_tithi_radians(   31, 3, 2004, 21.316666666666666, 5.5, &pdata);
//calculate_tithi_radians(   1 , 4, 2004, 21.716666666666665, 5.5, &pdata);
//calculate_tithi_radians(   2 , 4, 2004, 21.366666666666667, 5.5, &pdata);
//calculate_tithi_radians(   3 , 4, 2004, 20.333333333333332, 5.5, &pdata);
//calculate_tithi_radians(   4 , 4, 2004, 18.683333333333334, 5.5, &pdata);
//calculate_tithi_radians(   5 , 4, 2004, 16.533333333333335, 5.5, &pdata);
//calculate_tithi_radians(   6 , 4, 2004, 14.033333333333333, 5.5, &pdata);
//calculate_tithi_radians(   7 , 4, 2004, 11.3, 5.5, &pdata);
//calculate_tithi_radians(   8 , 4, 2004, 8.483333333333333, 5.5, &pdata);
//calculate_tithi_radians(   9 , 4, 2004, 5.666666666666667, 5.5, &pdata);
//calculate_tithi_radians(   10, 4, 2004, 2.966666666666667, 5.5, &pdata);
//calculate_tithi_radians(   11, 4, 2004, 0.4666666666666667, 5.5, &pdata);
//calculate_tithi_radians(   11, 4, 2004, 22.25, 5.5, &pdata);
//calculate_tithi_radians(   12, 4, 2004, 20.35, 5.5, &pdata);
//calculate_tithi_radians(   13, 4, 2004, 18.833333333333332, 5.5, &pdata);
//calculate_tithi_radians(   14, 4, 2004, 17.716666666666665, 5.5, &pdata);
//calculate_tithi_radians(   15, 4, 2004, 17.033333333333335, 5.5, &pdata);
//calculate_tithi_radians(   16, 4, 2004, 16.783333333333335, 5.5, &pdata);
//calculate_tithi_radians(   17, 4, 2004, 17.0, 5.5, &pdata);
//calculate_tithi_radians(   18, 4, 2004, 17.683333333333334, 5.5, &pdata);
//calculate_tithi_radians(   19, 4, 2004, 18.85, 5.5, &pdata);
//calculate_tithi_radians(   20, 4, 2004, 20.45, 5.5, &pdata);
//calculate_tithi_radians(   21, 4, 2004, 22.45, 5.5, &pdata);
//calculate_tithi_radians(   23, 4, 2004, 0.7833333333333333, 5.5, &pdata);
//calculate_tithi_radians(   24, 4, 2004, 3.3, 5.5, &pdata);
//calculate_tithi_radians(   25, 4, 2004, 5.85, 5.5, &pdata);
//calculate_tithi_radians(   26, 4, 2004, 8.233333333333333, 5.5, &pdata);
//calculate_tithi_radians(   27, 4, 2004, 10.233333333333333, 5.5, &pdata);
//calculate_tithi_radians(   28, 4, 2004, 11.683333333333334, 5.5, &pdata);
//calculate_tithi_radians(   29, 4, 2004, 12.433333333333334, 5.5, &pdata);
//calculate_tithi_radians(   30, 4, 2004, 12.433333333333334, 5.5, &pdata);
//calculate_tithi_radians(   1 , 5, 2004, 11.633333333333333, 5.5, &pdata);
//calculate_tithi_radians(   2 , 5, 2004, 10.1, 5.5, &pdata);
//calculate_tithi_radians(   3 , 5, 2004, 7.9, 5.5, &pdata);
//calculate_tithi_radians(   4 , 5, 2004, 5.183333333333334, 5.5, &pdata);
//calculate_tithi_radians(   5 , 5, 2004, 2.05, 5.5, &pdata);
//calculate_tithi_radians(   5 , 5, 2004, 22.666666666666668, 5.5, &pdata);
//calculate_tithi_radians(   6 , 5, 2004, 19.2, 5.5, &pdata);
//calculate_tithi_radians(   7 , 5, 2004, 15.766666666666667, 5.5, &pdata);
//calculate_tithi_radians(   8 , 5, 2004, 12.55, 5.5, &pdata);
//calculate_tithi_radians(   9 , 5, 2004, 9.65, 5.5, &pdata);
//calculate_tithi_radians(   10, 5, 2004, 7.2, 5.5, &pdata);
//calculate_tithi_radians(   11, 5, 2004, 5.3, 5.5, &pdata);
//calculate_tithi_radians(   12, 5, 2004, 3.9833333333333334, 5.5, &pdata);
//calculate_tithi_radians(   13, 5, 2004, 3.283333333333333, 5.5, &pdata);
//calculate_tithi_radians(   14, 5, 2004, 3.216666666666667, 5.5, &pdata);
//calculate_tithi_radians(   15, 5, 2004, 3.7333333333333334, 5.5, &pdata);
//calculate_tithi_radians(   16, 5, 2004, 4.783333333333333, 5.5, &pdata);
//calculate_tithi_radians(   17, 5, 2004, 6.283333333333333, 5.5, &pdata);
//calculate_tithi_radians(   18, 5, 2004, 8.166666666666666, 5.5, &pdata);

//calculate_tithi_radians( 14 , 10, 2023, 23.4, 5.5, &pdata);
//calculate_tithi_radians( 16 , 10, 2023, 0.5333333333333333, 5.5, &pdata);
//calculate_tithi_radians( 17 , 10, 2023, 1.2166666666666668, 5.5, &pdata);
//calculate_tithi_radians( 18 , 10, 2023, 1.4333333333333333, 5.5, &pdata);
//calculate_tithi_radians( 19 , 10, 2023, 1.2, 5.5, &pdata);
//calculate_tithi_radians( 20 , 10, 2023, 0.5166666666666667, 5.5, &pdata);
//calculate_tithi_radians( 20 , 10, 2023, 23.416666666666668, 5.5, &pdata);
//calculate_tithi_radians( 21 , 10, 2023, 21.883333333333333, 5.5, &pdata);
//calculate_tithi_radians( 22 , 10, 2023, 19.983333333333334, 5.5, &pdata);
//calculate_tithi_radians( 23 , 10, 2023, 17.75, 5.5, &pdata);
//calculate_tithi_radians( 24 , 10, 2023, 15.233333333333333, 5.5, &pdata);
//calculate_tithi_radians( 25 , 10, 2023, 12.533333333333333, 5.5, &pdata);
//calculate_tithi_radians( 26 , 10, 2023, 9.733333333333333, 5.5, &pdata);
//calculate_tithi_radians( 27 , 10, 2023, 6.95, 5.5, &pdata);
//calculate_tithi_radians( 28 , 10, 2023, 4.283333333333333, 5.5, &pdata);
//calculate_tithi_radians( 29 , 10, 2023, 1.8833333333333333, 5.5, &pdata);
//calculate_tithi_radians( 29 , 10, 2023, 23.883333333333333, 5.5, &pdata);
//calculate_tithi_radians( 30 , 10, 2023, 22.383333333333333, 5.5, &pdata);
//calculate_tithi_radians( 31 , 10, 2023, 21.5, 5.5, &pdata);
//calculate_tithi_radians( 1 , 11, 2023, 21.316666666666666, 5.5, &pdata);
//calculate_tithi_radians( 2 , 11, 2023, 21.866666666666667, 5.5, &pdata);
//calculate_tithi_radians( 3 , 11, 2023, 23.116666666666667, 5.5, &pdata);
//calculate_tithi_radians( 5 , 11, 2023, 0.9833333333333333, 5.5, &pdata);
//calculate_tithi_radians( 6 , 11, 2023, 3.3, 5.5, &pdata);
//calculate_tithi_radians( 7 , 11, 2023, 5.85, 5.5, &pdata);
//calculate_tithi_radians( 8 , 11, 2023, 8.383333333333333, 5.5, &pdata);
//calculate_tithi_radians( 9 , 11, 2023, 10.683333333333334, 5.5, &pdata);
//calculate_tithi_radians( 10 , 11, 2023, 12.583333333333334, 5.5, &pdata);
//calculate_tithi_radians( 11 , 11, 2023, 13.966666666666667, 5.5, &pdata);
//calculate_tithi_radians( 12 , 11, 2023, 14.75, 5.5, &pdata);
//calculate_tithi_radians( 13 , 11, 2023, 14.95, 5.5, &pdata);
//calculate_tithi_radians( 14 , 11, 2023, 14.6, 5.5, &pdata);
//calculate_tithi_radians( 15 , 11, 2023, 13.783333333333333, 5.5, &pdata);
//calculate_tithi_radians( 16 , 11, 2023, 12.566666666666666, 5.5, &pdata);
//calculate_tithi_radians( 17 , 11, 2023, 11.05, 5.5, &pdata);
//calculate_tithi_radians( 18 , 11, 2023, 9.3, 5.5, &pdata);
//calculate_tithi_radians( 19 , 11, 2023, 7.383333333333334, 5.5, &pdata);
//calculate_tithi_radians( 20 , 11, 2023, 5.35, 5.5, &pdata);
//calculate_tithi_radians( 21 , 11, 2023, 3.2666666666666666, 5.5, &pdata);
//calculate_tithi_radians( 22 , 11, 2023, 1.1666666666666667, 5.5, &pdata);
//calculate_tithi_radians( 22 , 11, 2023, 23.066666666666666, 5.5, &pdata);
//calculate_tithi_radians( 23 , 11, 2023, 21.033333333333335, 5.5, &pdata);
//calculate_tithi_radians( 24 , 11, 2023, 19.1, 5.5, &pdata);
//calculate_tithi_radians( 25 , 11, 2023, 17.366666666666667, 5.5, &pdata);
//calculate_tithi_radians( 26 , 11, 2023, 15.883333333333333, 5.5, &pdata);
//calculate_tithi_radians( 27 , 11, 2023, 14.766666666666667, 5.5, &pdata);
//calculate_tithi_radians( 28 , 11, 2023, 14.083333333333334, 5.5, &pdata);
//calculate_tithi_radians( 29 , 11, 2023, 13.95, 5.5, &pdata);
//calculate_tithi_radians( 30 , 11, 2023, 14.416666666666666, 5.5, &pdata);
//calculate_tithi_radians( 1 , 12, 2023, 15.516666666666667, 5.5, &pdata);
//calculate_tithi_radians( 2 , 12, 2023, 17.233333333333334, 5.5, &pdata);
//calculate_tithi_radians( 3 , 12, 2023, 19.45, 5.5, &pdata);
//calculate_tithi_radians( 4 , 12, 2023, 21.983333333333334, 5.5, &pdata);
//calculate_tithi_radians( 6 , 12, 2023, 0.6166666666666667, 5.5, &pdata);
//calculate_tithi_radians( 7 , 12, 2023, 3.066666666666667, 5.5, &pdata);
//calculate_tithi_radians( 8 , 12, 2023, 5.1, 5.5, &pdata);
//calculate_tithi_radians( 9 , 12, 2023, 6.516666666666667, 5.5, &pdata);
//calculate_tithi_radians( 10 , 12, 2023, 7.216666666666667, 5.5, &pdata);
//calculate_tithi_radians( 11 , 12, 2023, 7.166666666666667, 5.5, &pdata);
//calculate_tithi_radians( 12 , 12, 2023, 6.4, 5.5, &pdata);
//calculate_tithi_radians( 13 , 12, 2023, 5.016666666666667, 5.5, &pdata);
//calculate_tithi_radians( 14 , 12, 2023, 3.15, 5.5, &pdata);
//calculate_tithi_radians( 15 , 12, 2023, 0.9333333333333333, 5.5, &pdata);
//calculate_tithi_radians( 15 , 12, 2023, 22.5, 5.5, &pdata);
//calculate_tithi_radians( 16 , 12, 2023, 20.0, 5.5, &pdata);
//calculate_tithi_radians( 17 , 12, 2023, 17.55, 5.5, &pdata);
//calculate_tithi_radians( 18 , 12, 2023, 15.233333333333333, 5.5, &pdata);
//calculate_tithi_radians( 19 , 12, 2023, 13.116666666666667, 5.5, &pdata);
//calculate_tithi_radians( 20 , 12, 2023, 11.233333333333333, 5.5, &pdata);
//calculate_tithi_radians( 21 , 12, 2023, 9.616666666666667, 5.5, &pdata);
//calculate_tithi_radians( 22 , 12, 2023, 8.266666666666667, 5.5, &pdata);
//calculate_tithi_radians( 23 , 12, 2023, 7.2, 5.5, &pdata);
//calculate_tithi_radians( 24 , 12, 2023, 6.4, 5.5, &pdata);
//calculate_tithi_radians( 25 , 12, 2023, 5.916666666666667, 5.5, &pdata);
//calculate_tithi_radians( 26 , 12, 2023, 5.766666666666667, 5.5, &pdata);
//calculate_tithi_radians( 27 , 12, 2023, 6.033333333333333, 5.5, &pdata);
//calculate_tithi_radians( 28 , 12, 2023, 6.766666666666667, 5.5, &pdata);
//calculate_tithi_radians( 29 , 12, 2023, 7.983333333333333, 5.5, &pdata);
//calculate_tithi_radians( 30 , 12, 2023, 9.733333333333333, 5.5, &pdata);
//calculate_tithi_radians( 31 , 12, 2023, 11.916666666666666, 5.5, &pdata);
//calculate_tithi_radians( 1 , 1, 2024, 14.466666666666667, 5.5, &pdata);
//calculate_tithi_radians( 2 , 1, 2024, 17.166666666666668, 5.5, &pdata);
//calculate_tithi_radians( 3 , 1, 2024, 19.8, 5.5, &pdata);
//calculate_tithi_radians( 4 , 1, 2024, 22.066666666666666, 5.5, &pdata);
//calculate_tithi_radians( 5 , 1, 2024, 23.766666666666666, 5.5, &pdata);
//calculate_tithi_radians( 7 , 1, 2024, 0.6833333333333333, 5.5, &pdata);
//calculate_tithi_radians( 8 , 1, 2024, 0.7666666666666667, 5.5, &pdata);
//calculate_tithi_radians( 8 , 1, 2024, 23.983333333333334, 5.5, &pdata);
//calculate_tithi_radians( 9 , 1, 2024, 22.416666666666668, 5.5, &pdata);
//calculate_tithi_radians( 10 , 1, 2024, 20.183333333333334, 5.5, &pdata);
//calculate_tithi_radians( 11 , 1, 2024, 17.45, 5.5, &pdata);
//calculate_tithi_radians( 12 , 1, 2024, 14.383333333333333, 5.5, &pdata);
//calculate_tithi_radians( 13 , 1, 2024, 11.183333333333334, 5.5, &pdata);
//calculate_tithi_radians( 14 , 1, 2024, 8.0, 5.5, &pdata);
//calculate_tithi_radians( 15 , 1, 2024, 4.983333333333333, 5.5, &pdata);
//calculate_tithi_radians( 16 , 1, 2024, 2.2666666666666666, 5.5, &pdata);
//calculate_tithi_radians( 16 , 1, 2024, 23.966666666666665, 5.5, &pdata);
//calculate_tithi_radians( 17 , 1, 2024, 22.1, 5.5, &pdata);
//calculate_tithi_radians( 18 , 1, 2024, 20.75, 5.5, &pdata);
//calculate_tithi_radians( 19 , 1, 2024, 19.866666666666667, 5.5, &pdata);
//calculate_tithi_radians( 20 , 1, 2024, 19.433333333333334, 5.5, &pdata);
//calculate_tithi_radians( 21 , 1, 2024, 19.45, 5.5, &pdata);
//calculate_tithi_radians( 22 , 1, 2024, 19.85, 5.5, &pdata);
//calculate_tithi_radians( 23 , 1, 2024, 20.65, 5.5, &pdata);
//calculate_tithi_radians( 24 , 1, 2024, 21.833333333333332, 5.5, &pdata);
//calculate_tithi_radians( 25 , 1, 2024, 23.383333333333333, 5.5, &pdata);
//calculate_tithi_radians( 27 , 1, 2024, 1.3166666666666667, 5.5, &pdata);
//calculate_tithi_radians( 28 , 1, 2024, 3.6166666666666667, 5.5, &pdata);
//calculate_tithi_radians( 29 , 1, 2024, 6.183333333333334, 5.5, &pdata);
//calculate_tithi_radians( 30 , 1, 2024, 8.9, 5.5, &pdata);
//calculate_tithi_radians( 31 , 1, 2024, 11.6, 5.5, &pdata);
//calculate_tithi_radians( 1 , 2, 2024, 14.05, 5.5, &pdata);
//calculate_tithi_radians( 2 , 2, 2024, 16.033333333333335, 5.5, &pdata);
//calculate_tithi_radians( 3 , 2, 2024, 17.35, 5.5, &pdata);
//calculate_tithi_radians( 4 , 2, 2024, 17.816666666666666, 5.5, &pdata);
//calculate_tithi_radians( 5 , 2, 2024, 17.416666666666668, 5.5, &pdata);
//calculate_tithi_radians( 6 , 2, 2024, 16.116666666666667, 5.5, &pdata);
//calculate_tithi_radians( 7 , 2, 2024, 14.033333333333333, 5.5, &pdata);
//calculate_tithi_radians( 8 , 2, 2024, 11.283333333333333, 5.5, &pdata);
//calculate_tithi_radians( 9 , 2, 2024, 8.033333333333333, 5.5, &pdata);
//calculate_tithi_radians( 10 , 2, 2024, 4.466666666666667, 5.5, &pdata);
//calculate_tithi_radians( 11 , 2, 2024, 0.7833333333333333, 5.5, &pdata);
//calculate_tithi_radians( 11 , 2, 2024, 21.15, 5.5, &pdata);
//calculate_tithi_radians( 12 , 2, 2024, 17.733333333333334, 5.5, &pdata);
//calculate_tithi_radians( 13 , 2, 2024, 14.7, 5.5, &pdata);
//calculate_tithi_radians( 14 , 2, 2024, 12.15, 5.5, &pdata);
//calculate_tithi_radians( 15 , 2, 2024, 10.216666666666667, 5.5, &pdata);
//calculate_tithi_radians( 16 , 2, 2024, 8.9, 5.5, &pdata);
//calculate_tithi_radians( 17 , 2, 2024, 8.266666666666667, 5.5, &pdata);
//calculate_tithi_radians( 18 , 2, 2024, 8.25, 5.5, &pdata);
//calculate_tithi_radians( 19 , 2, 2024, 8.833333333333334, 5.5, &pdata);
//calculate_tithi_radians( 20 , 2, 2024, 9.916666666666666, 5.5, &pdata);
//calculate_tithi_radians( 21 , 2, 2024, 11.45, 5.5, &pdata);
//calculate_tithi_radians( 22 , 2, 2024, 13.35, 5.5, &pdata);
//calculate_tithi_radians( 23 , 2, 2024, 15.55, 5.5, &pdata);
//calculate_tithi_radians( 24 , 2, 2024, 18.0, 5.5, &pdata);
//calculate_tithi_radians( 25 , 2, 2024, 20.6, 5.5, &pdata);
//calculate_tithi_radians( 26 , 2, 2024, 23.266666666666666, 5.5, &pdata);
//calculate_tithi_radians( 28 , 2, 2024, 1.8833333333333333, 5.5, &pdata);
//calculate_tithi_radians( 29 , 2, 2024, 4.3, 5.5, &pdata);
//calculate_tithi_radians( 1 , 3, 2024, 6.366666666666666, 5.5, &pdata);
//calculate_tithi_radians( 2 , 3, 2024, 7.883333333333333, 5.5, &pdata);
//calculate_tithi_radians( 3 , 3, 2024, 8.75, 5.5, &pdata);
//calculate_tithi_radians( 4 , 3, 2024, 8.816666666666666, 5.5, &pdata);
//calculate_tithi_radians( 5 , 3, 2024, 8.066666666666666, 5.5, &pdata);
//calculate_tithi_radians( 6 , 3, 2024, 6.5, 5.5, &pdata);
//calculate_tithi_radians( 7 , 3, 2024, 4.216666666666667, 5.5, &pdata);
//calculate_tithi_radians( 8 , 3, 2024, 1.3166666666666667, 5.5, &pdata);
//calculate_tithi_radians( 8 , 3, 2024, 21.966666666666665, 5.5, &pdata);
//calculate_tithi_radians( 9 , 3, 2024, 18.283333333333335, 5.5, &pdata);
//calculate_tithi_radians( 10 , 3, 2024, 14.5, 5.5, &pdata);
//calculate_tithi_radians( 11 , 3, 2024, 10.75, 5.5, &pdata);
//calculate_tithi_radians( 12 , 3, 2024, 7.216666666666667, 5.5, &pdata);
//calculate_tithi_radians( 13 , 3, 2024, 4.066666666666666, 5.5, &pdata);
//calculate_tithi_radians( 14 , 3, 2024, 1.4333333333333333, 5.5, &pdata);
//calculate_tithi_radians( 14 , 3, 2024, 23.433333333333334, 5.5, &pdata);
//calculate_tithi_radians( 15 , 3, 2024, 22.15, 5.5, &pdata);
//calculate_tithi_radians( 16 , 3, 2024, 21.633333333333333, 5.5, &pdata);
//calculate_tithi_radians( 17 , 3, 2024, 21.883333333333333, 5.5, &pdata);
//calculate_tithi_radians( 18 , 3, 2024, 22.816666666666666, 5.5, &pdata);
//calculate_tithi_radians( 20 , 3, 2024, 0.36666666666666664, 5.5, &pdata);
//calculate_tithi_radians( 21 , 3, 2024, 2.3833333333333333, 5.5, &pdata);
//calculate_tithi_radians( 22 , 3, 2024, 4.733333333333333, 5.5, &pdata);
//calculate_tithi_radians( 23 , 3, 2024, 7.283333333333333, 5.5, &pdata);
//calculate_tithi_radians( 24 , 3, 2024, 9.916666666666666, 5.5, &pdata);
//calculate_tithi_radians( 25 , 3, 2024, 12.5, 5.5, &pdata);
//calculate_tithi_radians( 26 , 3, 2024, 14.916666666666666, 5.5, &pdata);
//calculate_tithi_radians( 27 , 3, 2024, 17.1, 5.5, &pdata);
//calculate_tithi_radians( 28 , 3, 2024, 18.933333333333334, 5.5, &pdata);
//calculate_tithi_radians( 29 , 3, 2024, 20.35, 5.5, &pdata);
//calculate_tithi_radians( 30 , 3, 2024, 21.233333333333334, 5.5, &pdata);
//calculate_tithi_radians( 31 , 3, 2024, 21.516666666666666, 5.5, &pdata);
//calculate_tithi_radians( 1 , 4, 2024, 21.15, 5.5, &pdata);
//calculate_tithi_radians( 2 , 4, 2024, 20.15, 5.5, &pdata);
//calculate_tithi_radians( 3 , 4, 2024, 18.483333333333334, 5.5, &pdata);
//calculate_tithi_radians( 4 , 4, 2024, 16.233333333333334, 5.5, &pdata);
//calculate_tithi_radians( 5 , 4, 2024, 13.466666666666667, 5.5, &pdata);
//calculate_tithi_radians( 6 , 4, 2024, 10.316666666666666, 5.5, &pdata);
//calculate_tithi_radians( 7 , 4, 2024, 6.9, 5.5, &pdata);
//calculate_tithi_radians( 8 , 4, 2024, 3.35, 5.5, &pdata);
//calculate_tithi_radians( 8 , 4, 2024, 23.833333333333332, 5.5, &pdata);
//calculate_tithi_radians( 9 , 4, 2024, 20.516666666666666, 5.5, &pdata);
//calculate_tithi_radians( 10 , 4, 2024, 17.533333333333335, 5.5, &pdata);
//calculate_tithi_radians( 11 , 4, 2024, 15.05, 5.5, &pdata);
//calculate_tithi_radians( 12 , 4, 2024, 13.2, 5.5, &pdata);
//calculate_tithi_radians( 13 , 4, 2024, 12.066666666666666, 5.5, &pdata);
//calculate_tithi_radians( 14 , 4, 2024, 11.733333333333333, 5.5, &pdata);
//calculate_tithi_radians( 15 , 4, 2024, 12.183333333333334, 5.5, &pdata);
//calculate_tithi_radians( 16 , 4, 2024, 13.4, 5.5, &pdata);
//calculate_tithi_radians( 17 , 4, 2024, 15.233333333333333, 5.5, &pdata);
//calculate_tithi_radians( 18 , 4, 2024, 17.516666666666666, 5.5, &pdata);
//calculate_tithi_radians( 19 , 4, 2024, 20.083333333333332, 5.5, &pdata);
//calculate_tithi_radians( 20 , 4, 2024, 22.683333333333334, 5.5, &pdata);
//calculate_tithi_radians( 22 , 4, 2024, 1.1833333333333333, 5.5, &pdata);
//calculate_tithi_radians( 23 , 4, 2024, 3.4166666666666665, 5.5, &pdata);
//calculate_tithi_radians( 24 , 4, 2024, 5.3, 5.5, &pdata);
//calculate_tithi_radians( 25 , 4, 2024, 6.766666666666667, 5.5, &pdata);
//calculate_tithi_radians( 26 , 4, 2024, 7.766666666666667, 5.5, &pdata);
//calculate_tithi_radians( 27 , 4, 2024, 8.3, 5.5, &pdata);
//calculate_tithi_radians( 28 , 4, 2024, 8.35, 5.5, &pdata);
//calculate_tithi_radians( 29 , 4, 2024, 7.95, 5.5, &pdata);
//calculate_tithi_radians( 30 , 4, 2024, 7.083333333333333, 5.5, &pdata);
//calculate_tithi_radians( 1 , 5, 2024, 5.766666666666667, 5.5, &pdata);
//calculate_tithi_radians( 2 , 5, 2024, 4.016666666666667, 5.5, &pdata);
//calculate_tithi_radians( 3 , 5, 2024, 1.8833333333333333, 5.5, &pdata);
//calculate_tithi_radians( 3 , 5, 2024, 23.4, 5.5, &pdata);
//calculate_tithi_radians( 4 , 5, 2024, 20.633333333333333, 5.5, &pdata);
//calculate_tithi_radians( 5 , 5, 2024, 17.7, 5.5, &pdata);
//calculate_tithi_radians( 6 , 5, 2024, 14.666666666666666, 5.5, &pdata);
//calculate_tithi_radians( 7 , 5, 2024, 11.666666666666666, 5.5, &pdata);
//calculate_tithi_radians( 8 , 5, 2024, 8.85, 5.5, &pdata);
//calculate_tithi_radians( 9 , 5, 2024, 6.35, 5.5, &pdata);
//calculate_tithi_radians( 10 , 5, 2024, 4.3, 5.5, &pdata);
//calculate_tithi_radians( 11 , 5, 2024, 2.8333333333333335, 5.5, &pdata);
//calculate_tithi_radians( 12 , 5, 2024, 2.066666666666667, 5.5, &pdata);
//calculate_tithi_radians( 13 , 5, 2024, 2.066666666666667, 5.5, &pdata);
//calculate_tithi_radians( 14 , 5, 2024, 2.8333333333333335, 5.5, &pdata);
//calculate_tithi_radians( 15 , 5, 2024, 4.316666666666666, 5.5, &pdata);
//calculate_tithi_radians( 16 , 5, 2024, 6.383333333333334, 5.5, &pdata);
//calculate_tithi_radians( 17 , 5, 2024, 8.816666666666666, 5.5, &pdata);
//calculate_tithi_radians( 18 , 5, 2024, 11.383333333333333, 5.5, &pdata);
//calculate_tithi_radians( 19 , 5, 2024, 13.833333333333334, 5.5, &pdata);
//calculate_tithi_radians( 20 , 5, 2024, 15.966666666666667, 5.5, &pdata);
//calculate_tithi_radians( 21 , 5, 2024, 17.65, 5.5, &pdata);
//calculate_tithi_radians( 22 , 5, 2024, 18.8, 5.5, &pdata);
//calculate_tithi_radians( 23 , 5, 2024, 19.366666666666667, 5.5, &pdata);
//calculate_tithi_radians( 24 , 5, 2024, 19.416666666666668, 5.5, &pdata);
//calculate_tithi_radians( 25 , 5, 2024, 18.966666666666665, 5.5, &pdata);
//calculate_tithi_radians( 26 , 5, 2024, 18.1, 5.5, &pdata);
//calculate_tithi_radians( 27 , 5, 2024, 16.883333333333333, 5.5, &pdata);
//calculate_tithi_radians( 28 , 5, 2024, 15.383333333333333, 5.5, &pdata);
//calculate_tithi_radians( 29 , 5, 2024, 13.65, 5.5, &pdata);
//calculate_tithi_radians( 30 , 5, 2024, 11.733333333333333, 5.5, &pdata);
//calculate_tithi_radians( 31 , 5, 2024, 9.633333333333333, 5.5, &pdata);
//calculate_tithi_radians( 1 , 6, 2024, 7.4, 5.5, &pdata);
//calculate_tithi_radians( 2 , 6, 2024, 5.066666666666666, 5.5, &pdata);
//calculate_tithi_radians( 3 , 6, 2024, 2.6833333333333336, 5.5, &pdata);
//calculate_tithi_radians( 4 , 6, 2024, 0.3, 5.5, &pdata);
//calculate_tithi_radians( 4 , 6, 2024, 22.016666666666666, 5.5, &pdata);
//calculate_tithi_radians( 5 , 6, 2024, 19.916666666666668, 5.5, &pdata);
//calculate_tithi_radians( 6 , 6, 2024, 18.116666666666667, 5.5, &pdata);
//calculate_tithi_radians( 7 , 6, 2024, 16.75, 5.5, &pdata);
//calculate_tithi_radians( 8 , 6, 2024, 15.916666666666666, 5.5, &pdata);
//calculate_tithi_radians( 9 , 6, 2024, 15.733333333333333, 5.5, &pdata);
//calculate_tithi_radians( 10 , 6, 2024, 16.25, 5.5, &pdata);
//calculate_tithi_radians( 11 , 6, 2024, 17.45, 5.5, &pdata);
//calculate_tithi_radians( 12 , 6, 2024, 19.266666666666666, 5.5, &pdata);
//calculate_tithi_radians( 13 , 6, 2024, 21.55, 5.5, &pdata);
//calculate_tithi_radians( 15 , 6, 2024, 0.06666666666666667, 5.5, &pdata);
//calculate_tithi_radians( 16 , 6, 2024, 2.533333333333333, 5.5, &pdata);
//calculate_tithi_radians( 17 , 6, 2024, 4.716666666666667, 5.5, &pdata);
//calculate_tithi_radians( 18 , 6, 2024, 6.416666666666667, 5.5, &pdata);
//calculate_tithi_radians( 19 , 6, 2024, 7.466666666666667, 5.5, &pdata);
//calculate_tithi_radians( 20 , 6, 2024, 7.833333333333333, 5.5, &pdata);
//calculate_tithi_radians( 21 , 6, 2024, 7.516666666666667, 5.5, &pdata);
//calculate_tithi_radians( 22 , 6, 2024, 6.616666666666667, 5.5, &pdata);
//calculate_tithi_radians( 23 , 6, 2024, 5.216666666666667, 5.5, &pdata);
//calculate_tithi_radians( 24 , 6, 2024, 3.4333333333333336, 5.5, &pdata);
//calculate_tithi_radians( 25 , 6, 2024, 1.3833333333333333, 5.5, &pdata);
//calculate_tithi_radians( 25 , 6, 2024, 23.183333333333334, 5.5, &pdata);
//calculate_tithi_radians( 26 , 6, 2024, 20.916666666666668, 5.5, &pdata);
//calculate_tithi_radians( 27 , 6, 2024, 18.65, 5.5, &pdata);
//calculate_tithi_radians( 28 , 6, 2024, 16.45, 5.5, &pdata);
//calculate_tithi_radians( 29 , 6, 2024, 14.333333333333334, 5.5, &pdata);
//calculate_tithi_radians( 30 , 6, 2024, 12.316666666666666, 5.5, &pdata);
//calculate_tithi_radians( 1 , 7, 2024, 10.433333333333334, 5.5, &pdata);
//calculate_tithi_radians( 2 , 7, 2024, 8.7, 5.5, &pdata);
//calculate_tithi_radians( 3 , 7, 2024, 7.166666666666667, 5.5, &pdata);
//calculate_tithi_radians( 4 , 7, 2024, 5.9, 5.5, &pdata);
//calculate_tithi_radians( 5 , 7, 2024, 4.966666666666667, 5.5, &pdata);
//calculate_tithi_radians( 6 , 7, 2024, 4.45, 5.5, &pdata);
//calculate_tithi_radians( 7 , 7, 2024, 4.433333333333334, 5.5, &pdata);
//calculate_tithi_radians( 8 , 7, 2024, 4.983333333333333, 5.5, &pdata);
//calculate_tithi_radians( 9 , 7, 2024, 6.133333333333334, 5.5, &pdata);
//calculate_tithi_radians( 10 , 7, 2024, 7.866666666666667, 5.5, &pdata);
//calculate_tithi_radians( 11 , 7, 2024, 10.05, 5.5, &pdata);
//calculate_tithi_radians( 12 , 7, 2024, 12.55, 5.5, &pdata);
//calculate_tithi_radians( 13 , 7, 2024, 15.083333333333334, 5.5, &pdata);
//calculate_tithi_radians( 14 , 7, 2024, 17.433333333333334, 5.5, &pdata);
//calculate_tithi_radians( 15 , 7, 2024, 19.316666666666666, 5.5, &pdata);
//calculate_tithi_radians( 16 , 7, 2024, 20.566666666666666, 5.5, &pdata);
//calculate_tithi_radians( 17 , 7, 2024, 21.033333333333335, 5.5, &pdata);
//calculate_tithi_radians( 18 , 7, 2024, 20.733333333333334, 5.5, &pdata);
//calculate_tithi_radians( 19 , 7, 2024, 19.683333333333334, 5.5, &pdata);
//calculate_tithi_radians( 20 , 7, 2024, 17.983333333333334, 5.5, &pdata);
//calculate_tithi_radians( 21 , 7, 2024, 15.766666666666667, 5.5, &pdata);
//calculate_tithi_radians( 22 , 7, 2024, 13.183333333333334, 5.5, &pdata);
//calculate_tithi_radians( 23 , 7, 2024, 10.383333333333333, 5.5, &pdata);
//calculate_tithi_radians( 24 , 7, 2024, 7.5, 5.5, &pdata);
//calculate_tithi_radians( 25 , 7, 2024, 4.666666666666667, 5.5, &pdata);
//calculate_tithi_radians( 26 , 7, 2024, 1.9666666666666668, 5.5, &pdata);
//calculate_tithi_radians( 26 , 7, 2024, 23.5, 5.5, &pdata);
//calculate_tithi_radians( 27 , 7, 2024, 21.316666666666666, 5.5, &pdata);
//calculate_tithi_radians( 28 , 7, 2024, 19.45, 5.5, &pdata);
//calculate_tithi_radians( 29 , 7, 2024, 17.933333333333334, 5.5, &pdata);
//calculate_tithi_radians( 30 , 7, 2024, 16.75, 5.5, &pdata);
//calculate_tithi_radians( 31 , 7, 2024, 15.916666666666666, 5.5, &pdata);
//calculate_tithi_radians( 1 , 8, 2024, 15.483333333333333, 5.5, &pdata);
//calculate_tithi_radians( 2 , 8, 2024, 15.433333333333334, 5.5, &pdata);
//calculate_tithi_radians( 3 , 8, 2024, 15.833333333333334, 5.5, &pdata);
//calculate_tithi_radians( 4 , 8, 2024, 16.7, 5.5, &pdata);
//calculate_tithi_radians( 5 , 8, 2024, 18.05, 5.5, &pdata);
//calculate_tithi_radians( 6 , 8, 2024, 19.866666666666667, 5.5, &pdata);
//calculate_tithi_radians( 7 , 8, 2024, 22.1, 5.5, &pdata);
//calculate_tithi_radians( 9 , 8, 2024, 0.6, 5.5, &pdata);
//calculate_tithi_radians( 10 , 8, 2024, 3.2333333333333334, 5.5, &pdata);
//calculate_tithi_radians( 11 , 8, 2024, 5.75, 5.5, &pdata);
//calculate_tithi_radians( 12 , 8, 2024, 7.916666666666667, 5.5, &pdata);
//calculate_tithi_radians( 13 , 8, 2024, 9.516666666666667, 5.5, &pdata);
//calculate_tithi_radians( 14 , 8, 2024, 10.383333333333333, 5.5, &pdata);
//calculate_tithi_radians( 15 , 8, 2024, 10.433333333333334, 5.5, &pdata);
//calculate_tithi_radians( 16 , 8, 2024, 9.65, 5.5, &pdata);
//calculate_tithi_radians( 17 , 8, 2024, 8.083333333333334, 5.5, &pdata);
//calculate_tithi_radians( 18 , 8, 2024, 5.85, 5.5, &pdata);
//calculate_tithi_radians( 19 , 8, 2024, 3.066666666666667, 5.5, &pdata);
//calculate_tithi_radians( 19 , 8, 2024, 23.916666666666668, 5.5, &pdata);
//calculate_tithi_radians( 20 , 8, 2024, 20.55, 5.5, &pdata);
//calculate_tithi_radians( 21 , 8, 2024, 17.116666666666667, 5.5, &pdata);
//calculate_tithi_radians( 22 , 8, 2024, 13.766666666666667, 5.5, &pdata);
//calculate_tithi_radians( 23 , 8, 2024, 10.65, 5.5, &pdata);
//calculate_tithi_radians( 24 , 8, 2024, 7.866666666666667, 5.5, &pdata);
//calculate_tithi_radians( 25 , 8, 2024, 5.5, 5.5, &pdata);
//calculate_tithi_radians( 26 , 8, 2024, 3.65, 5.5, &pdata);
//calculate_tithi_radians( 27 , 8, 2024, 2.3166666666666664, 5.5, &pdata);
//calculate_tithi_radians( 28 , 8, 2024, 1.55, 5.5, &pdata);
//calculate_tithi_radians( 29 , 8, 2024, 1.3166666666666667, 5.5, &pdata);
//calculate_tithi_radians( 30 , 8, 2024, 1.6166666666666667, 5.5, &pdata);
//calculate_tithi_radians( 31 , 8, 2024, 2.4166666666666665, 5.5, &pdata);
//calculate_tithi_radians( 1 , 9, 2024, 3.6833333333333336, 5.5, &pdata);
//calculate_tithi_radians( 2 , 9, 2024, 5.35, 5.5, &pdata);
//calculate_tithi_radians( 3 , 9, 2024, 7.416666666666667, 5.5, &pdata);
//calculate_tithi_radians( 4 , 9, 2024, 9.783333333333333, 5.5, &pdata);
//calculate_tithi_radians( 5 , 9, 2024, 12.35, 5.5, &pdata);
//calculate_tithi_radians( 6 , 9, 2024, 15.016666666666667, 5.5, &pdata);
//calculate_tithi_radians( 7 , 9, 2024, 17.616666666666667, 5.5, &pdata);
//calculate_tithi_radians( 8 , 9, 2024, 19.966666666666665, 5.5, &pdata);
//calculate_tithi_radians( 9 , 9, 2024, 21.883333333333333, 5.5, &pdata);
//calculate_tithi_radians( 10 , 9, 2024, 23.2, 5.5, &pdata);
//calculate_tithi_radians( 11 , 9, 2024, 23.766666666666666, 5.5, &pdata);
//calculate_tithi_radians( 12 , 9, 2024, 23.55, 5.5, &pdata);
//calculate_tithi_radians( 13 , 9, 2024, 22.5, 5.5, &pdata);
//calculate_tithi_radians( 14 , 9, 2024, 20.683333333333334, 5.5, &pdata);
//calculate_tithi_radians( 15 , 9, 2024, 18.2, 5.5, &pdata);
//calculate_tithi_radians( 16 , 9, 2024, 15.166666666666666, 5.5, &pdata);
//calculate_tithi_radians( 17 , 9, 2024, 11.733333333333333, 5.5, &pdata);
//calculate_tithi_radians( 18 , 9, 2024, 8.066666666666666, 5.5, &pdata);
//calculate_tithi_radians( 19 , 9, 2024, 4.316666666666666, 5.5, &pdata);
//calculate_tithi_radians( 20 , 9, 2024, 0.6666666666666666, 5.5, &pdata);
//calculate_tithi_radians( 20 , 9, 2024, 21.25, 5.5, &pdata);
//calculate_tithi_radians( 21 , 9, 2024, 18.233333333333334, 5.5, &pdata);
//calculate_tithi_radians( 22 , 9, 2024, 15.716666666666667, 5.5, &pdata);
//calculate_tithi_radians( 23 , 9, 2024, 13.833333333333334, 5.5, &pdata);
//calculate_tithi_radians( 24 , 9, 2024, 12.65, 5.5, &pdata);
//calculate_tithi_radians( 25 , 9, 2024, 12.183333333333334, 5.5, &pdata);
//calculate_tithi_radians( 26 , 9, 2024, 12.416666666666666, 5.5, &pdata);
//calculate_tithi_radians( 27 , 9, 2024, 13.333333333333334, 5.5, &pdata);
//calculate_tithi_radians( 28 , 9, 2024, 14.833333333333334, 5.5, &pdata);
//calculate_tithi_radians( 29 , 9, 2024, 16.783333333333335, 5.5, &pdata);
//calculate_tithi_radians( 30 , 9, 2024, 19.1, 5.5, &pdata);
//calculate_tithi_radians( 1 , 10, 2024, 21.65, 5.5, &pdata);
//calculate_tithi_radians( 3 , 10, 2024, 0.31666666666666665, 5.5, &pdata);
//calculate_tithi_radians( 4 , 10, 2024, 2.966666666666667, 5.5, &pdata);
//calculate_tithi_radians( 5 , 10, 2024, 5.5, 5.5, &pdata);
//calculate_tithi_radians( 6 , 10, 2024, 7.816666666666666, 5.5, &pdata);
//calculate_tithi_radians( 7 , 10, 2024, 9.783333333333333, 5.5, &pdata);
//calculate_tithi_radians( 8 , 10, 2024, 11.3, 5.5, &pdata);
//calculate_tithi_radians( 9 , 10, 2024, 12.233333333333333, 5.5, &pdata);
//calculate_tithi_radians( 10 , 10, 2024, 12.516666666666667, 5.5, &pdata);
//calculate_tithi_radians( 11 , 10, 2024, 12.1, 5.5, &pdata);
//calculate_tithi_radians( 12 , 10, 2024, 10.966666666666667, 5.5, &pdata);
//calculate_tithi_radians( 13 , 10, 2024, 9.133333333333333, 5.5, &pdata);
//calculate_tithi_radians( 14 , 10, 2024, 6.683333333333334, 5.5, &pdata);
//calculate_tithi_radians( 15 , 10, 2024, 3.7, 5.5, &pdata);
//calculate_tithi_radians( 16 , 10, 2024, 0.31666666666666665, 5.5, &pdata);
//calculate_tithi_radians( 16 , 10, 2024, 20.666666666666668, 5.5, &pdata);
//calculate_tithi_radians( 17 , 10, 2024, 16.933333333333334, 5.5, &pdata);
//calculate_tithi_radians( 18 , 10, 2024, 13.25, 5.5, &pdata);
//calculate_tithi_radians( 19 , 10, 2024, 9.8, 5.5, &pdata);
//calculate_tithi_radians( 20 , 10, 2024, 6.766666666666667, 5.5, &pdata);
//calculate_tithi_radians( 21 , 10, 2024, 4.283333333333333, 5.5, &pdata);
//calculate_tithi_radians( 22 , 10, 2024, 2.4833333333333334, 5.5, &pdata);
//calculate_tithi_radians( 23 , 10, 2024, 1.4833333333333334, 5.5, &pdata);
//calculate_tithi_radians( 24 , 10, 2024, 1.3166666666666667, 5.5, &pdata);
//calculate_tithi_radians( 25 , 10, 2024, 1.9666666666666668, 5.5, &pdata);
//calculate_tithi_radians( 26 , 10, 2024, 3.3833333333333333, 5.5, &pdata);
//calculate_tithi_radians( 27 , 10, 2024, 5.4, 5.5, &pdata);
//calculate_tithi_radians( 28 , 10, 2024, 7.833333333333333, 5.5, &pdata);
//calculate_tithi_radians( 29 , 10, 2024, 10.516666666666667, 5.5, &pdata);
//calculate_tithi_radians( 30 , 10, 2024, 13.25, 5.5, &pdata);
//calculate_tithi_radians( 31 , 10, 2024, 15.883333333333333, 5.5, &pdata);
//calculate_tithi_radians( 1 , 11, 2024, 18.266666666666666, 5.5, &pdata);
//calculate_tithi_radians( 2 , 11, 2024, 20.366666666666667, 5.5, &pdata);
//calculate_tithi_radians( 3 , 11, 2024, 22.083333333333332, 5.5, &pdata);
//calculate_tithi_radians( 4 , 11, 2024, 23.4, 5.5, &pdata);
//calculate_tithi_radians( 6 , 11, 2024, 0.2833333333333333, 5.5, &pdata);
//calculate_tithi_radians( 7 , 11, 2024, 0.6833333333333333, 5.5, &pdata);
//calculate_tithi_radians( 8 , 11, 2024, 0.5833333333333334, 5.5, &pdata);
//calculate_tithi_radians( 8 , 11, 2024, 23.933333333333334, 5.5, &pdata);
//calculate_tithi_radians( 9 , 11, 2024, 22.75, 5.5, &pdata);
//calculate_tithi_radians( 10 , 11, 2024, 21.016666666666666, 5.5, &pdata);
//calculate_tithi_radians( 11 , 11, 2024, 18.766666666666666, 5.5, &pdata);
//calculate_tithi_radians( 12 , 11, 2024, 16.083333333333332, 5.5, &pdata);
//calculate_tithi_radians( 13 , 11, 2024, 13.016666666666667, 5.5, &pdata);
//calculate_tithi_radians( 14 , 11, 2024, 9.716666666666667, 5.5, &pdata);
//calculate_tithi_radians( 15 , 11, 2024, 6.316666666666666, 5.5, &pdata);
//calculate_tithi_radians( 16 , 11, 2024, 2.966666666666667, 5.5, &pdata);
//calculate_tithi_radians( 16 , 11, 2024, 23.833333333333332, 5.5, &pdata);
//calculate_tithi_radians( 17 , 11, 2024, 21.1, 5.5, &pdata);
//calculate_tithi_radians( 18 , 11, 2024, 18.933333333333334, 5.5, &pdata);
//calculate_tithi_radians( 19 , 11, 2024, 17.466666666666665, 5.5, &pdata);
//calculate_tithi_radians( 20 , 11, 2024, 16.816666666666666, 5.5, &pdata);
//calculate_tithi_radians( 21 , 11, 2024, 17.05, 5.5, &pdata);
//calculate_tithi_radians( 22 , 11, 2024, 18.133333333333333, 5.5, &pdata);
//calculate_tithi_radians( 23 , 11, 2024, 19.95, 5.5, &pdata);
//calculate_tithi_radians( 24 , 11, 2024, 22.333333333333332, 5.5, &pdata);
//calculate_tithi_radians( 26 , 11, 2024, 1.0166666666666666, 5.5, &pdata);
//calculate_tithi_radians( 27 , 11, 2024, 3.783333333333333, 5.5, &pdata);
//calculate_tithi_radians( 28 , 11, 2024, 6.383333333333334, 5.5, &pdata);
//calculate_tithi_radians( 29 , 11, 2024, 8.666666666666666, 5.5, &pdata);
//calculate_tithi_radians( 30 , 11, 2024, 10.5, 5.5, &pdata);
//calculate_tithi_radians( 1 , 12, 2024, 11.85, 5.5, &pdata);
//calculate_tithi_radians( 2 , 12, 2024, 12.716666666666667, 5.5, &pdata);
//calculate_tithi_radians( 3 , 12, 2024, 13.15, 5.5, &pdata);
//calculate_tithi_radians( 4 , 12, 2024, 13.166666666666666, 5.5, &pdata);
//calculate_tithi_radians( 5 , 12, 2024, 12.816666666666666, 5.5, &pdata);
//calculate_tithi_radians( 6 , 12, 2024, 12.133333333333333, 5.5, &pdata);
//calculate_tithi_radians( 7 , 12, 2024, 11.1, 5.5, &pdata);
//calculate_tithi_radians( 8 , 12, 2024, 9.733333333333333, 5.5, &pdata);
//calculate_tithi_radians( 9 , 12, 2024, 8.05, 5.5, &pdata);
//calculate_tithi_radians( 10 , 12, 2024, 6.016666666666667, 5.5, &pdata);
//calculate_tithi_radians( 11 , 12, 2024, 3.716666666666667, 5.5, &pdata);
//calculate_tithi_radians( 12 , 12, 2024, 1.15, 5.5, &pdata);
//calculate_tithi_radians( 12 , 12, 2024, 22.433333333333334, 5.5, &pdata);
//calculate_tithi_radians( 13 , 12, 2024, 19.666666666666668, 5.5, &pdata);
//calculate_tithi_radians( 14 , 12, 2024, 16.966666666666665, 5.5, &pdata);
//calculate_tithi_radians( 15 , 12, 2024, 14.516666666666667, 5.5, &pdata);
//calculate_tithi_radians( 16 , 12, 2024, 12.45, 5.5, &pdata);
//calculate_tithi_radians( 17 , 12, 2024, 10.933333333333334, 5.5, &pdata);
//calculate_tithi_radians( 18 , 12, 2024, 10.1, 5.5, &pdata);
//calculate_tithi_radians( 19 , 12, 2024, 10.05, 5.5, &pdata);
//calculate_tithi_radians( 20 , 12, 2024, 10.816666666666666, 5.5, &pdata);
//calculate_tithi_radians( 21 , 12, 2024, 12.35, 5.5, &pdata);
//calculate_tithi_radians( 22 , 12, 2024, 14.533333333333333, 5.5, &pdata);
//calculate_tithi_radians( 23 , 12, 2024, 17.116666666666667, 5.5, &pdata);
//calculate_tithi_radians( 24 , 12, 2024, 19.866666666666667, 5.5, &pdata);
//calculate_tithi_radians( 25 , 12, 2024, 22.483333333333334, 5.5, &pdata);
//calculate_tithi_radians( 27 , 12, 2024, 0.7333333333333333, 5.5, &pdata);
//calculate_tithi_radians( 28 , 12, 2024, 2.4333333333333336, 5.5, &pdata);
//calculate_tithi_radians( 29 , 12, 2024, 3.533333333333333, 5.5, &pdata);
//calculate_tithi_radians( 30 , 12, 2024, 4.016666666666667, 5.5, &pdata);
//calculate_tithi_radians( 31 , 12, 2024, 3.9333333333333336, 5.5, &pdata);
//calculate_tithi_radians( 1 , 1, 2025, 3.3666666666666667, 5.5, &pdata);
//calculate_tithi_radians( 2 , 1, 2025, 2.4, 5.5, &pdata);
//calculate_tithi_radians( 3 , 1, 2025, 1.1333333333333333, 5.5, &pdata);
//calculate_tithi_radians( 3 , 1, 2025, 23.65, 5.5, &pdata);
//calculate_tithi_radians( 4 , 1, 2025, 22.016666666666666, 5.5, &pdata);
//calculate_tithi_radians( 5 , 1, 2025, 20.25, 5.5, &pdata);
//calculate_tithi_radians( 6 , 1, 2025, 18.383333333333333, 5.5, &pdata);
//calculate_tithi_radians( 7 , 1, 2025, 16.45, 5.5, &pdata);
//calculate_tithi_radians( 8 , 1, 2025, 14.433333333333334, 5.5, &pdata);
//calculate_tithi_radians( 9 , 1, 2025, 12.366666666666667, 5.5, &pdata);
//calculate_tithi_radians( 10 , 1, 2025, 10.316666666666666, 5.5, &pdata);
//calculate_tithi_radians( 11 , 1, 2025, 8.35, 5.5, &pdata);
//calculate_tithi_radians( 12 , 1, 2025, 6.55, 5.5, &pdata);
//calculate_tithi_radians( 13 , 1, 2025, 5.05, 5.5, &pdata);
//calculate_tithi_radians( 14 , 1, 2025, 3.9333333333333336, 5.5, &pdata);
//calculate_tithi_radians( 15 , 1, 2025, 3.35, 5.5, &pdata);
//calculate_tithi_radians( 16 , 1, 2025, 3.3833333333333333, 5.5, &pdata);
//calculate_tithi_radians( 17 , 1, 2025, 4.1, 5.5, &pdata);
//calculate_tithi_radians( 18 , 1, 2025, 5.5, 5.5, &pdata);
//calculate_tithi_radians( 19 , 1, 2025, 7.5, 5.5, &pdata);
//calculate_tithi_radians( 20 , 1, 2025, 9.966666666666667, 5.5, &pdata);
//calculate_tithi_radians( 21 , 1, 2025, 12.65, 5.5, &pdata);
//calculate_tithi_radians( 22 , 1, 2025, 15.3, 5.5, &pdata);
//calculate_tithi_radians( 23 , 1, 2025, 17.616666666666667, 5.5, &pdata);
//calculate_tithi_radians( 24 , 1, 2025, 19.416666666666668, 5.5, &pdata);
//calculate_tithi_radians( 25 , 1, 2025, 20.533333333333335, 5.5, &pdata);
//calculate_tithi_radians( 26 , 1, 2025, 20.916666666666668, 5.5, &pdata);
//calculate_tithi_radians( 27 , 1, 2025, 20.566666666666666, 5.5, &pdata);
//calculate_tithi_radians( 28 , 1, 2025, 19.6, 5.5, &pdata);
//calculate_tithi_radians( 29 , 1, 2025, 18.083333333333332, 5.5, &pdata);
//calculate_tithi_radians( 30 , 1, 2025, 16.166666666666668, 5.5, &pdata);
//calculate_tithi_radians( 31 , 1, 2025, 13.983333333333333, 5.5, &pdata);
//calculate_tithi_radians( 1 , 2, 2025, 11.633333333333333, 5.5, &pdata);
//calculate_tithi_radians( 2 , 2, 2025, 9.233333333333333, 5.5, &pdata);
//calculate_tithi_radians( 3 , 2, 2025, 6.866666666666667, 5.5, &pdata);
//calculate_tithi_radians( 4 , 2, 2025, 4.616666666666667, 5.5, &pdata);
//calculate_tithi_radians( 5 , 2, 2025, 2.5, 5.5, &pdata);
//calculate_tithi_radians( 6 , 2, 2025, 0.5833333333333334, 5.5, &pdata);
//calculate_tithi_radians( 6 , 2, 2025, 22.883333333333333, 5.5, &pdata);
//calculate_tithi_radians( 7 , 2, 2025, 21.433333333333334, 5.5, &pdata);
//calculate_tithi_radians( 8 , 2, 2025, 20.266666666666666, 5.5, &pdata);
//calculate_tithi_radians( 9 , 2, 2025, 19.416666666666668, 5.5, &pdata);
//calculate_tithi_radians( 10 , 2, 2025, 18.95, 5.5, &pdata);
//calculate_tithi_radians( 11 , 2, 2025, 18.916666666666668, 5.5, &pdata);
//calculate_tithi_radians( 12 , 2, 2025, 19.383333333333333, 5.5, &pdata);
//calculate_tithi_radians( 13 , 2, 2025, 20.35, 5.5, &pdata);
//calculate_tithi_radians( 14 , 2, 2025, 21.866666666666667, 5.5, &pdata);
//calculate_tithi_radians( 15 , 2, 2025, 23.866666666666667, 5.5, &pdata);
//calculate_tithi_radians( 17 , 2, 2025, 2.2666666666666666, 5.5, &pdata);
//calculate_tithi_radians( 18 , 2, 2025, 4.883333333333333, 5.5, &pdata);
//calculate_tithi_radians( 19 , 2, 2025, 7.533333333333333, 5.5, &pdata);
//calculate_tithi_radians( 20 , 2, 2025, 9.966666666666667, 5.5, &pdata);
//calculate_tithi_radians( 21 , 2, 2025, 11.95, 5.5, &pdata);
//calculate_tithi_radians( 22 , 2, 2025, 13.316666666666666, 5.5, &pdata);
//calculate_tithi_radians( 23 , 2, 2025, 13.933333333333334, 5.5, &pdata);
//calculate_tithi_radians( 24 , 2, 2025, 13.75, 5.5, &pdata);
//calculate_tithi_radians( 25 , 2, 2025, 12.783333333333333, 5.5, &pdata);
//calculate_tithi_radians( 26 , 2, 2025, 11.133333333333333, 5.5, &pdata);
//calculate_tithi_radians( 27 , 2, 2025, 8.9, 5.5, &pdata);
//calculate_tithi_radians( 28 , 2, 2025, 6.233333333333333, 5.5, &pdata);
//calculate_tithi_radians( 1 , 3, 2025, 3.2666666666666666, 5.5, &pdata);
//calculate_tithi_radians( 2 , 3, 2025, 0.15, 5.5, &pdata);
//calculate_tithi_radians( 2 , 3, 2025, 21.033333333333335, 5.5, &pdata);
//calculate_tithi_radians( 3 , 3, 2025, 18.033333333333335, 5.5, &pdata);
//calculate_tithi_radians( 4 , 3, 2025, 15.266666666666667, 5.5, &pdata);
//calculate_tithi_radians( 5 , 3, 2025, 12.85, 5.5, &pdata);
//calculate_tithi_radians( 6 , 3, 2025, 10.85, 5.5, &pdata);
//calculate_tithi_radians( 7 , 3, 2025, 9.3, 5.5, &pdata);
//calculate_tithi_radians( 8 , 3, 2025, 8.266666666666667, 5.5, &pdata);
//calculate_tithi_radians( 9 , 3, 2025, 7.75, 5.5, &pdata);
//calculate_tithi_radians( 10 , 3, 2025, 7.733333333333333, 5.5, &pdata);
//calculate_tithi_radians( 11 , 3, 2025, 8.233333333333333, 5.5, &pdata);
//calculate_tithi_radians( 12 , 3, 2025, 9.183333333333334, 5.5, &pdata);
//calculate_tithi_radians( 13 , 3, 2025, 10.6, 5.5, &pdata);
//calculate_tithi_radians( 14 , 3, 2025, 12.4, 5.5, &pdata);
//calculate_tithi_radians( 15 , 3, 2025, 14.55, 5.5, &pdata);
//calculate_tithi_radians( 16 , 3, 2025, 16.966666666666665, 5.5, &pdata);
//calculate_tithi_radians( 17 , 3, 2025, 19.55, 5.5, &pdata);
//calculate_tithi_radians( 18 , 3, 2025, 22.15, 5.5, &pdata);
//calculate_tithi_radians( 20 , 3, 2025, 0.6166666666666667, 5.5, &pdata);
//calculate_tithi_radians( 21 , 3, 2025, 2.75, 5.5, &pdata);
//calculate_tithi_radians( 22 , 3, 2025, 4.4, 5.5, &pdata);
//calculate_tithi_radians( 23 , 3, 2025, 5.383333333333334, 5.5, &pdata);
//calculate_tithi_radians( 24 , 3, 2025, 5.633333333333333, 5.5, &pdata);
//calculate_tithi_radians( 25 , 3, 2025, 5.083333333333333, 5.5, &pdata);
//calculate_tithi_radians( 26 , 3, 2025, 3.75, 5.5, &pdata);
//calculate_tithi_radians( 27 , 3, 2025, 1.7166666666666668, 5.5, &pdata);
//calculate_tithi_radians( 27 , 3, 2025, 23.05, 5.5, &pdata);
//calculate_tithi_radians( 28 , 3, 2025, 19.916666666666668, 5.5, &pdata);
//calculate_tithi_radians( 29 , 3, 2025, 16.45, 5.5, &pdata);
//calculate_tithi_radians( 30 , 3, 2025, 12.816666666666666, 5.5, &pdata);
//calculate_tithi_radians( 31 , 3, 2025, 9.183333333333334, 5.5, &pdata);
//calculate_tithi_radians( 1 , 4, 2025, 5.7, 5.5, &pdata);
//calculate_tithi_radians( 2 , 4, 2025, 2.533333333333333, 5.5, &pdata);
//calculate_tithi_radians( 2 , 4, 2025, 23.816666666666666, 5.5, &pdata);
//calculate_tithi_radians( 3 , 4, 2025, 21.683333333333334, 5.5, &pdata);
//calculate_tithi_radians( 4 , 4, 2025, 20.2, 5.5, &pdata);
//calculate_tithi_radians( 5 , 4, 2025, 19.433333333333334, 5.5, &pdata);
//calculate_tithi_radians( 6 , 4, 2025, 19.383333333333333, 5.5, &pdata);
//calculate_tithi_radians( 7 , 4, 2025, 20.0, 5.5, &pdata);
//calculate_tithi_radians( 8 , 4, 2025, 21.216666666666665, 5.5, &pdata);
//calculate_tithi_radians( 9 , 4, 2025, 22.916666666666668, 5.5, &pdata);
//calculate_tithi_radians( 11 , 4, 2025, 1.0, 5.5, &pdata);
//calculate_tithi_radians( 12 , 4, 2025, 3.35, 5.5, &pdata);
//calculate_tithi_radians( 13 , 4, 2025, 5.85, 5.5, &pdata);
//calculate_tithi_radians( 14 , 4, 2025, 8.416666666666666, 5.5, &pdata);
//calculate_tithi_radians( 15 , 4, 2025, 10.916666666666666, 5.5, &pdata);
//calculate_tithi_radians( 16 , 4, 2025, 13.283333333333333, 5.5, &pdata);
//calculate_tithi_radians( 17 , 4, 2025, 15.383333333333333, 5.5, &pdata);
//calculate_tithi_radians( 18 , 4, 2025, 17.116666666666667, 5.5, &pdata);
//calculate_tithi_radians( 19 , 4, 2025, 18.366666666666667, 5.5, &pdata);
//calculate_tithi_radians( 20 , 4, 2025, 19.0, 5.5, &pdata);
//calculate_tithi_radians( 21 , 4, 2025, 18.966666666666665, 5.5, &pdata);
//calculate_tithi_radians( 22 , 4, 2025, 18.216666666666665, 5.5, &pdata);
//calculate_tithi_radians( 23 , 4, 2025, 16.716666666666665, 5.5, &pdata);
//calculate_tithi_radians( 24 , 4, 2025, 14.533333333333333, 5.5, &pdata);
//calculate_tithi_radians( 25 , 4, 2025, 11.733333333333333, 5.5, &pdata);
//calculate_tithi_radians( 26 , 4, 2025, 8.45, 5.5, &pdata);
//calculate_tithi_radians( 27 , 4, 2025, 4.816666666666666, 5.5, &pdata);
//calculate_tithi_radians( 28 , 4, 2025, 1.0, 5.5, &pdata);
//calculate_tithi_radians( 28 , 4, 2025, 21.183333333333334, 5.5, &pdata);
//calculate_tithi_radians( 29 , 4, 2025, 17.516666666666666, 5.5, &pdata);
//calculate_tithi_radians( 30 , 4, 2025, 14.2, 5.5, &pdata);
//calculate_tithi_radians( 1 , 5, 2025, 11.383333333333333, 5.5, &pdata);
//calculate_tithi_radians( 2 , 5, 2025, 9.233333333333333, 5.5, &pdata);
//calculate_tithi_radians( 3 , 5, 2025, 7.866666666666667, 5.5, &pdata);
//calculate_tithi_radians( 4 , 5, 2025, 7.316666666666666, 5.5, &pdata);
//calculate_tithi_radians( 5 , 5, 2025, 7.583333333333333, 5.5, &pdata);
//calculate_tithi_radians( 6 , 5, 2025, 8.633333333333333, 5.5, &pdata);
//calculate_tithi_radians( 7 , 5, 2025, 10.316666666666666, 5.5, &pdata);
//calculate_tithi_radians( 8 , 5, 2025, 12.483333333333333, 5.5, &pdata);
//calculate_tithi_radians( 9 , 5, 2025, 14.933333333333334, 5.5, &pdata);
//calculate_tithi_radians( 10 , 5, 2025, 17.5, 5.5, &pdata);
//calculate_tithi_radians( 11 , 5, 2025, 20.033333333333335, 5.5, &pdata);
//calculate_tithi_radians( 12 , 5, 2025, 22.416666666666668, 5.5, &pdata);
//calculate_tithi_radians( 14 , 5, 2025, 0.5833333333333334, 5.5, &pdata);
//calculate_tithi_radians( 15 , 5, 2025, 2.4833333333333334, 5.5, &pdata);
//calculate_tithi_radians( 16 , 5, 2025, 4.05, 5.5, &pdata);
//calculate_tithi_radians( 17 , 5, 2025, 5.216666666666667, 5.5, &pdata);
//calculate_tithi_radians( 18 , 5, 2025, 5.95, 5.5, &pdata);
//calculate_tithi_radians( 19 , 5, 2025, 6.183333333333334, 5.5, &pdata);
//calculate_tithi_radians( 20 , 5, 2025, 5.85, 5.5, &pdata);
//calculate_tithi_radians( 21 , 5, 2025, 4.916666666666667, 5.5, &pdata);
//calculate_tithi_radians( 22 , 5, 2025, 3.3666666666666667, 5.5, &pdata);
//calculate_tithi_radians( 23 , 5, 2025, 1.2, 5.5, &pdata);
//calculate_tithi_radians( 23 , 5, 2025, 22.483333333333334, 5.5, &pdata);
//calculate_tithi_radians( 24 , 5, 2025, 19.333333333333332, 5.5, &pdata);
//calculate_tithi_radians( 25 , 5, 2025, 15.85, 5.5, &pdata);
//calculate_tithi_radians( 26 , 5, 2025, 12.183333333333334, 5.5, &pdata);
//calculate_tithi_radians( 27 , 5, 2025, 8.533333333333333, 5.5, &pdata);
//calculate_tithi_radians( 28 , 5, 2025, 5.033333333333333, 5.5, &pdata);
//calculate_tithi_radians( 29 , 5, 2025, 1.9, 5.5, &pdata);
//calculate_tithi_radians( 29 , 5, 2025, 23.3, 5.5, &pdata);
//calculate_tithi_radians( 30 , 5, 2025, 21.366666666666667, 5.5, &pdata);
//calculate_tithi_radians( 31 , 5, 2025, 20.25, 5.5, &pdata);
//calculate_tithi_radians( 1 , 6, 2025, 19.983333333333334, 5.5, &pdata);
//calculate_tithi_radians( 2 , 6, 2025, 20.583333333333332, 5.5, &pdata);
//calculate_tithi_radians( 3 , 6, 2025, 21.933333333333334, 5.5, &pdata);
//calculate_tithi_radians( 4 , 6, 2025, 23.9, 5.5, &pdata);
//calculate_tithi_radians( 6 , 6, 2025, 2.2666666666666666, 5.5, &pdata);
//calculate_tithi_radians( 7 , 6, 2025, 4.8, 5.5, &pdata);
//calculate_tithi_radians( 8 , 6, 2025, 7.3, 5.5, &pdata);
//calculate_tithi_radians( 9 , 6, 2025, 9.6, 5.5, &pdata);
//calculate_tithi_radians( 10 , 6, 2025, 11.583333333333334, 5.5, &pdata);
//calculate_tithi_radians( 11 , 6, 2025, 13.216666666666667, 5.5, &pdata);
//calculate_tithi_radians( 12 , 6, 2025, 14.466666666666667, 5.5, &pdata);
//calculate_tithi_radians( 13 , 6, 2025, 15.316666666666666, 5.5, &pdata);
//calculate_tithi_radians( 14 , 6, 2025, 15.783333333333333, 5.5, &pdata);
//calculate_tithi_radians( 15 , 6, 2025, 15.85, 5.5, &pdata);
//calculate_tithi_radians( 16 , 6, 2025, 15.516666666666667, 5.5, &pdata);
//calculate_tithi_radians( 17 , 6, 2025, 14.766666666666667, 5.5, &pdata);
//calculate_tithi_radians( 18 , 6, 2025, 13.583333333333334, 5.5, &pdata);
//calculate_tithi_radians( 19 , 6, 2025, 11.916666666666666, 5.5, &pdata);
//calculate_tithi_radians( 20 , 6, 2025, 9.816666666666666, 5.5, &pdata);
//calculate_tithi_radians( 21 , 6, 2025, 7.3, 5.5, &pdata);
//calculate_tithi_radians( 22 , 6, 2025, 4.45, 5.5, &pdata);
//calculate_tithi_radians( 23 , 6, 2025, 1.3666666666666667, 5.5, &pdata);
//calculate_tithi_radians( 23 , 6, 2025, 22.15, 5.5, &pdata);
//calculate_tithi_radians( 24 , 6, 2025, 18.983333333333334, 5.5, &pdata);
//calculate_tithi_radians( 25 , 6, 2025, 16.016666666666666, 5.5, &pdata);
//calculate_tithi_radians( 26 , 6, 2025, 13.4, 5.5, &pdata);
//calculate_tithi_radians( 27 , 6, 2025, 11.316666666666666, 5.5, &pdata);
//calculate_tithi_radians( 28 , 6, 2025, 9.9, 5.5, &pdata);
//calculate_tithi_radians( 29 , 6, 2025, 9.233333333333333, 5.5, &pdata);
//calculate_tithi_radians( 30 , 6, 2025, 9.4, 5.5, &pdata);
//calculate_tithi_radians( 1 , 7, 2025, 10.333333333333334, 5.5, &pdata);
//calculate_tithi_radians( 2 , 7, 2025, 11.966666666666667, 5.5, &pdata);
//calculate_tithi_radians( 3 , 7, 2025, 14.116666666666667, 5.5, &pdata);
//calculate_tithi_radians( 4 , 7, 2025, 16.533333333333335, 5.5, &pdata);
//calculate_tithi_radians( 5 , 7, 2025, 18.983333333333334, 5.5, &pdata);
//calculate_tithi_radians( 6 , 7, 2025, 21.25, 5.5, &pdata);
//calculate_tithi_radians( 7 , 7, 2025, 23.166666666666668, 5.5, &pdata);
//calculate_tithi_radians( 9 , 7, 2025, 0.6333333333333333, 5.5, &pdata);
//calculate_tithi_radians( 10 , 7, 2025, 1.6166666666666667, 5.5, &pdata);
//calculate_tithi_radians( 11 , 7, 2025, 2.1, 5.5, &pdata);
//calculate_tithi_radians( 12 , 7, 2025, 2.1333333333333333, 5.5, &pdata);
//calculate_tithi_radians( 13 , 7, 2025, 1.7666666666666666, 5.5, &pdata);
//calculate_tithi_radians( 14 , 7, 2025, 1.0333333333333334, 5.5, &pdata);
//calculate_tithi_radians( 14 , 7, 2025, 23.983333333333334, 5.5, &pdata);
//calculate_tithi_radians( 15 , 7, 2025, 22.65, 5.5, &pdata);
//calculate_tithi_radians( 16 , 7, 2025, 21.033333333333335, 5.5, &pdata);
//calculate_tithi_radians( 17 , 7, 2025, 19.15, 5.5, &pdata);
//calculate_tithi_radians( 18 , 7, 2025, 17.033333333333335, 5.5, &pdata);
//calculate_tithi_radians( 19 , 7, 2025, 14.7, 5.5, &pdata);
//calculate_tithi_radians( 20 , 7, 2025, 12.216666666666667, 5.5, &pdata);
//calculate_tithi_radians( 21 , 7, 2025, 9.65, 5.5, &pdata);
//calculate_tithi_radians( 22 , 7, 2025, 7.083333333333333, 5.5, &pdata);
//calculate_tithi_radians( 23 , 7, 2025, 4.65, 5.5, &pdata);
//calculate_tithi_radians( 24 , 7, 2025, 2.466666666666667, 5.5, &pdata);
//calculate_tithi_radians( 25 , 7, 2025, 0.6666666666666666, 5.5, &pdata);
//calculate_tithi_radians( 25 , 7, 2025, 23.383333333333333, 5.5, &pdata);
//calculate_tithi_radians( 26 , 7, 2025, 22.7, 5.5, &pdata);
//calculate_tithi_radians( 27 , 7, 2025, 22.7, 5.5, &pdata);
//calculate_tithi_radians( 28 , 7, 2025, 23.4, 5.5, &pdata);
//calculate_tithi_radians( 30 , 7, 2025, 0.7666666666666667, 5.5, &pdata);
//calculate_tithi_radians( 31 , 7, 2025, 2.6833333333333336, 5.5, &pdata);
//calculate_tithi_radians( 1 , 8, 2025, 4.966666666666667, 5.5, &pdata);
//calculate_tithi_radians( 2 , 8, 2025, 7.383333333333334, 5.5, &pdata);
//calculate_tithi_radians( 3 , 8, 2025, 9.7, 5.5, &pdata);
//calculate_tithi_radians( 4 , 8, 2025, 11.683333333333334, 5.5, &pdata);
//calculate_tithi_radians( 5 , 8, 2025, 13.2, 5.5, &pdata);
//calculate_tithi_radians( 6 , 8, 2025, 14.133333333333333, 5.5, &pdata);
//calculate_tithi_radians( 7 , 8, 2025, 14.466666666666667, 5.5, &pdata);
//calculate_tithi_radians( 8 , 8, 2025, 14.2, 5.5, &pdata);
//calculate_tithi_radians( 9 , 8, 2025, 13.4, 5.5, &pdata);
//calculate_tithi_radians( 10 , 8, 2025, 12.166666666666666, 5.5, &pdata);
//calculate_tithi_radians( 11 , 8, 2025, 10.55, 5.5, &pdata);
//calculate_tithi_radians( 12 , 8, 2025, 8.683333333333334, 5.5, &pdata);
//calculate_tithi_radians( 13 , 8, 2025, 6.6, 5.5, &pdata);
//calculate_tithi_radians( 14 , 8, 2025, 4.383333333333334, 5.5, &pdata);
//calculate_tithi_radians( 15 , 8, 2025, 2.1166666666666667, 5.5, &pdata);
//calculate_tithi_radians( 15 , 8, 2025, 23.833333333333332, 5.5, &pdata);
//calculate_tithi_radians( 16 , 8, 2025, 21.566666666666666, 5.5, &pdata);
//calculate_tithi_radians( 17 , 8, 2025, 19.4, 5.5, &pdata);
//calculate_tithi_radians( 18 , 8, 2025, 17.366666666666667, 5.5, &pdata);
//calculate_tithi_radians( 19 , 8, 2025, 15.533333333333333, 5.5, &pdata);
//calculate_tithi_radians( 20 , 8, 2025, 13.966666666666667, 5.5, &pdata);
//calculate_tithi_radians( 21 , 8, 2025, 12.75, 5.5, &pdata);
//calculate_tithi_radians( 22 , 8, 2025, 11.933333333333334, 5.5, &pdata);
//calculate_tithi_radians( 23 , 8, 2025, 11.6, 5.5, &pdata);
//calculate_tithi_radians( 24 , 8, 2025, 11.8, 5.5, &pdata);
//calculate_tithi_radians( 25 , 8, 2025, 12.583333333333334, 5.5, &pdata);
//calculate_tithi_radians( 26 , 8, 2025, 13.9, 5.5, &pdata);
//calculate_tithi_radians( 27 , 8, 2025, 15.733333333333333, 5.5, &pdata);
//calculate_tithi_radians( 28 , 8, 2025, 17.95, 5.5, &pdata);
//calculate_tithi_radians( 29 , 8, 2025, 20.366666666666667, 5.5, &pdata);
//calculate_tithi_radians( 30 , 8, 2025, 22.766666666666666, 5.5, &pdata);
//calculate_tithi_radians( 1 , 9, 2025, 0.95, 5.5, &pdata);
//calculate_tithi_radians( 2 , 9, 2025, 2.716666666666667, 5.5, &pdata);
//calculate_tithi_radians( 3 , 9, 2025, 3.8833333333333333, 5.5, &pdata);
//calculate_tithi_radians( 4 , 9, 2025, 4.366666666666666, 5.5, &pdata);
//calculate_tithi_radians( 5 , 9, 2025, 4.133333333333334, 5.5, &pdata);
//calculate_tithi_radians( 6 , 9, 2025, 3.216666666666667, 5.5, &pdata);
//calculate_tithi_radians( 7 , 9, 2025, 1.6833333333333333, 5.5, &pdata);
//calculate_tithi_radians( 7 , 9, 2025, 23.633333333333333, 5.5, &pdata);
//calculate_tithi_radians( 8 , 9, 2025, 21.2, 5.5, &pdata);
//calculate_tithi_radians( 9 , 9, 2025, 18.483333333333334, 5.5, &pdata);
//calculate_tithi_radians( 10 , 9, 2025, 15.633333333333333, 5.5, &pdata);
//calculate_tithi_radians( 11 , 9, 2025, 12.75, 5.5, &pdata);
//calculate_tithi_radians( 12 , 9, 2025, 9.966666666666667, 5.5, &pdata);
//calculate_tithi_radians( 13 , 9, 2025, 7.383333333333334, 5.5, &pdata);
//calculate_tithi_radians( 14 , 9, 2025, 5.066666666666666, 5.5, &pdata);
//calculate_tithi_radians( 15 , 9, 2025, 3.1, 5.5, &pdata);
//calculate_tithi_radians( 16 , 9, 2025, 1.5166666666666666, 5.5, &pdata);
//calculate_tithi_radians( 17 , 9, 2025, 0.36666666666666664, 5.5, &pdata);
//calculate_tithi_radians( 17 , 9, 2025, 23.65, 5.5, &pdata);
//calculate_tithi_radians( 18 , 9, 2025, 23.4, 5.5, &pdata);
//calculate_tithi_radians( 19 , 9, 2025, 23.616666666666667, 5.5, &pdata);
//calculate_tithi_radians( 21 , 9, 2025, 0.2833333333333333, 5.5, &pdata);
//calculate_tithi_radians( 22 , 9, 2025, 1.3833333333333333, 5.5, &pdata);
//calculate_tithi_radians( 23 , 9, 2025, 2.9333333333333336, 5.5, &pdata);
//calculate_tithi_radians( 24 , 9, 2025, 4.85, 5.5, &pdata);
//calculate_tithi_radians( 25 , 9, 2025, 7.1, 5.5, &pdata);
//calculate_tithi_radians( 26 , 9, 2025, 9.55, 5.5, &pdata);
//calculate_tithi_radians( 27 , 9, 2025, 12.05, 5.5, &pdata);
//calculate_tithi_radians( 28 , 9, 2025, 14.45, 5.5, &pdata);
//calculate_tithi_radians( 29 , 9, 2025, 16.516666666666666, 5.5, &pdata);
//calculate_tithi_radians( 30 , 9, 2025, 18.1, 5.5, &pdata);
//calculate_tithi_radians( 1 , 10, 2025, 19.016666666666666, 5.5, &pdata);
//calculate_tithi_radians( 2 , 10, 2025, 19.183333333333334, 5.5, &pdata);
//calculate_tithi_radians( 3 , 10, 2025, 18.55, 5.5, &pdata);
//calculate_tithi_radians( 4 , 10, 2025, 17.15, 5.5, &pdata);
//calculate_tithi_radians( 5 , 10, 2025, 15.066666666666666, 5.5, &pdata);
//calculate_tithi_radians( 6 , 10, 2025, 12.383333333333333, 5.5, &pdata);
//calculate_tithi_radians( 7 , 10, 2025, 9.283333333333333, 5.5, &pdata);
//calculate_tithi_radians( 8 , 10, 2025, 5.883333333333333, 5.5, &pdata);
//calculate_tithi_radians( 9 , 10, 2025, 2.3666666666666667, 5.5, &pdata);
//calculate_tithi_radians( 9 , 10, 2025, 22.9, 5.5, &pdata);
//calculate_tithi_radians( 10 , 10, 2025, 19.633333333333333, 5.5, &pdata);
//calculate_tithi_radians( 11 , 10, 2025, 16.716666666666665, 5.5, &pdata);
//calculate_tithi_radians( 12 , 10, 2025, 14.283333333333333, 5.5, &pdata);
//calculate_tithi_radians( 13 , 10, 2025, 12.4, 5.5, &pdata);
//calculate_tithi_radians( 14 , 10, 2025, 11.15, 5.5, &pdata);
//calculate_tithi_radians( 15 , 10, 2025, 10.55, 5.5, &pdata);
//calculate_tithi_radians( 16 , 10, 2025, 10.583333333333334, 5.5, &pdata);
//calculate_tithi_radians( 17 , 10, 2025, 11.2, 5.5, &pdata);
//calculate_tithi_radians( 18 , 10, 2025, 12.316666666666666, 5.5, &pdata);
//calculate_tithi_radians( 19 , 10, 2025, 13.85, 5.5, &pdata);
//calculate_tithi_radians( 20 , 10, 2025, 15.75, 5.5, &pdata);
//calculate_tithi_radians( 21 , 10, 2025, 17.9, 5.5, &pdata);
//calculate_tithi_radians( 22 , 10, 2025, 20.266666666666666, 5.5, &pdata);
//calculate_tithi_radians( 23 , 10, 2025, 22.766666666666666, 5.5, &pdata);
//calculate_tithi_radians( 25 , 10, 2025, 1.3166666666666667, 5.5, &pdata);
//calculate_tithi_radians( 26 , 10, 2025, 3.8, 5.5, &pdata);
//calculate_tithi_radians( 27 , 10, 2025, 6.083333333333333, 5.5, &pdata);
//calculate_tithi_radians( 28 , 10, 2025, 7.983333333333333, 5.5, &pdata);
//calculate_tithi_radians( 29 , 10, 2025, 9.383333333333333, 5.5, &pdata);
//calculate_tithi_radians( 30 , 10, 2025, 10.1, 5.5, &pdata);
//calculate_tithi_radians( 31 , 10, 2025, 10.05, 5.5, &pdata);
//calculate_tithi_radians( 1 , 11, 2025, 9.183333333333334, 5.5, &pdata);
//calculate_tithi_radians( 2 , 11, 2025, 7.516666666666667, 5.5, &pdata);
//calculate_tithi_radians( 3 , 11, 2025, 5.116666666666666, 5.5, &pdata);
//calculate_tithi_radians( 4 , 11, 2025, 2.1, 5.5, &pdata);
//calculate_tithi_radians( 4 , 11, 2025, 22.6, 5.5, &pdata);
//calculate_tithi_radians( 5 , 11, 2025, 18.816666666666666, 5.5, &pdata);
//calculate_tithi_radians( 6 , 11, 2025, 14.9, 5.5, &pdata);
//calculate_tithi_radians( 7 , 11, 2025, 11.083333333333334, 5.5, &pdata);
//calculate_tithi_radians( 8 , 11, 2025, 7.533333333333333, 5.5, &pdata);
//calculate_tithi_radians( 9 , 11, 2025, 4.416666666666667, 5.5, &pdata);
//calculate_tithi_radians( 10 , 11, 2025, 1.9166666666666665, 5.5, &pdata);
//calculate_tithi_radians( 11 , 11, 2025, 0.13333333333333333, 5.5, &pdata);
//calculate_tithi_radians( 11 , 11, 2025, 23.15, 5.5, &pdata);
//calculate_tithi_radians( 12 , 11, 2025, 22.966666666666665, 5.5, &pdata);
//calculate_tithi_radians( 13 , 11, 2025, 23.566666666666666, 5.5, &pdata);
//calculate_tithi_radians( 15 , 11, 2025, 0.8166666666666667, 5.5, &pdata);
//calculate_tithi_radians( 16 , 11, 2025, 2.6166666666666667, 5.5, &pdata);
//calculate_tithi_radians( 17 , 11, 2025, 4.783333333333333, 5.5, &pdata);
//calculate_tithi_radians( 18 , 11, 2025, 7.2, 5.5, &pdata);
//calculate_tithi_radians( 19 , 11, 2025, 9.716666666666667, 5.5, &pdata);
//calculate_tithi_radians( 20 , 11, 2025, 12.266666666666667, 5.5, &pdata);
//calculate_tithi_radians( 21 , 11, 2025, 14.783333333333333, 5.5, &pdata);
//calculate_tithi_radians( 22 , 11, 2025, 17.183333333333334, 5.5, &pdata);
//calculate_tithi_radians( 23 , 11, 2025, 19.4, 5.5, &pdata);
//calculate_tithi_radians( 24 , 11, 2025, 21.366666666666667, 5.5, &pdata);
//calculate_tithi_radians( 25 , 11, 2025, 22.95, 5.5, &pdata);
//calculate_tithi_radians( 27 , 11, 2025, 0.03333333333333333, 5.5, &pdata);
//calculate_tithi_radians( 28 , 11, 2025, 0.5, 5.5, &pdata);
//calculate_tithi_radians( 29 , 11, 2025, 0.25, 5.5, &pdata);
//calculate_tithi_radians( 29 , 11, 2025, 23.25, 5.5, &pdata);
//calculate_tithi_radians( 30 , 11, 2025, 21.483333333333334, 5.5, &pdata);
//calculate_tithi_radians( 1 , 12, 2025, 19.016666666666666, 5.5, &pdata);
//calculate_tithi_radians( 2 , 12, 2025, 15.95, 5.5, &pdata);
//calculate_tithi_radians( 3 , 12, 2025, 12.433333333333334, 5.5, &pdata);
//calculate_tithi_radians( 4 , 12, 2025, 8.616666666666667, 5.5, &pdata);
//calculate_tithi_radians( 5 , 12, 2025, 4.716666666666667, 5.5, &pdata);
//calculate_tithi_radians( 6 , 12, 2025, 0.9166666666666666, 5.5, &pdata);
//calculate_tithi_radians( 6 , 12, 2025, 21.416666666666668, 5.5, &pdata);
//calculate_tithi_radians( 7 , 12, 2025, 18.416666666666668, 5.5, &pdata);
//calculate_tithi_radians( 8 , 12, 2025, 16.05, 5.5, &pdata);
//calculate_tithi_radians( 9 , 12, 2025, 14.483333333333333, 5.5, &pdata);
//calculate_tithi_radians( 10 , 12, 2025, 13.766666666666667, 5.5, &pdata);
//calculate_tithi_radians( 11 , 12, 2025, 13.95, 5.5, &pdata);
//calculate_tithi_radians( 12 , 12, 2025, 14.95, 5.5, &pdata);
//calculate_tithi_radians( 13 , 12, 2025, 16.633333333333333, 5.5, &pdata);
//calculate_tithi_radians( 14 , 12, 2025, 18.816666666666666, 5.5, &pdata);
//calculate_tithi_radians( 15 , 12, 2025, 21.333333333333332, 5.5, &pdata);
//calculate_tithi_radians( 16 , 12, 2025, 23.95, 5.5, &pdata);
//calculate_tithi_radians( 18 , 12, 2025, 2.533333333333333, 5.5, &pdata);
//calculate_tithi_radians( 19 , 12, 2025, 4.983333333333333, 5.5, &pdata);
//calculate_tithi_radians( 20 , 12, 2025, 7.216666666666667, 5.5, &pdata);
//calculate_tithi_radians( 21 , 12, 2025, 9.183333333333334, 5.5, &pdata);
//calculate_tithi_radians( 22 , 12, 2025, 10.85, 5.5, &pdata);
//calculate_tithi_radians( 23 , 12, 2025, 12.216666666666667, 5.5, &pdata);
//calculate_tithi_radians( 24 , 12, 2025, 13.183333333333334, 5.5, &pdata);
//calculate_tithi_radians( 25 , 12, 2025, 13.7, 5.5, &pdata);
//calculate_tithi_radians( 26 , 12, 2025, 13.716666666666667, 5.5, &pdata);
//calculate_tithi_radians( 27 , 12, 2025, 13.166666666666666, 5.5, &pdata);
//calculate_tithi_radians( 28 , 12, 2025, 11.983333333333333, 5.5, &pdata);
//calculate_tithi_radians( 29 , 12, 2025, 10.2, 5.5, &pdata);
//calculate_tithi_radians( 30 , 12, 2025, 7.85, 5.5, &pdata);
//calculate_tithi_radians( 31 , 12, 2025, 5.0, 5.5, &pdata);
//calculate_tithi_radians( 1 , 1, 2026, 1.8, 5.5, &pdata);
//calculate_tithi_radians( 1 , 1, 2026, 22.366666666666667, 5.5, &pdata);
//calculate_tithi_radians( 2 , 1, 2026, 18.883333333333333, 5.5, &pdata);
//calculate_tithi_radians( 3 , 1, 2026, 15.533333333333333, 5.5, &pdata);
//calculate_tithi_radians( 4 , 1, 2026, 12.5, 5.5, &pdata);
//calculate_tithi_radians( 5 , 1, 2026, 9.933333333333334, 5.5, &pdata);
//calculate_tithi_radians( 6 , 1, 2026, 8.016666666666667, 5.5, &pdata);
//calculate_tithi_radians( 7 , 1, 2026, 6.866666666666667, 5.5, &pdata);
//calculate_tithi_radians( 8 , 1, 2026, 6.55, 5.5, &pdata);
//calculate_tithi_radians( 9 , 1, 2026, 7.083333333333333, 5.5, &pdata);
//calculate_tithi_radians( 10 , 1, 2026, 8.383333333333333, 5.5, &pdata);
//calculate_tithi_radians( 11 , 1, 2026, 10.333333333333334, 5.5, &pdata);
//calculate_tithi_radians( 12 , 1, 2026, 12.7, 5.5, &pdata);
//calculate_tithi_radians( 13 , 1, 2026, 15.3, 5.5, &pdata);
//calculate_tithi_radians( 14 , 1, 2026, 17.866666666666667, 5.5, &pdata);
//calculate_tithi_radians( 15 , 1, 2026, 20.266666666666666, 5.5, &pdata);
//calculate_tithi_radians( 16 , 1, 2026, 22.35, 5.5, &pdata);
//calculate_tithi_radians( 18 , 1, 2026, 0.06666666666666667, 5.5, &pdata);
//calculate_tithi_radians( 19 , 1, 2026, 1.35, 5.5, &pdata);
//calculate_tithi_radians( 20 , 1, 2026, 2.2333333333333334, 5.5, &pdata);
//calculate_tithi_radians( 21 , 1, 2026, 2.7, 5.5, &pdata);
//calculate_tithi_radians( 22 , 1, 2026, 2.783333333333333, 5.5, &pdata);
//calculate_tithi_radians( 23 , 1, 2026, 2.466666666666667, 5.5, &pdata);
//calculate_tithi_radians( 24 , 1, 2026, 1.7666666666666666, 5.5, &pdata);
//calculate_tithi_radians( 25 , 1, 2026, 0.6666666666666666, 5.5, &pdata);
//calculate_tithi_radians( 25 , 1, 2026, 23.166666666666668, 5.5, &pdata);
//calculate_tithi_radians( 26 , 1, 2026, 21.3, 5.5, &pdata);
//calculate_tithi_radians( 27 , 1, 2026, 19.083333333333332, 5.5, &pdata);
//calculate_tithi_radians( 28 , 1, 2026, 16.6, 5.5, &pdata);
//calculate_tithi_radians( 29 , 1, 2026, 13.916666666666666, 5.5, &pdata);
//calculate_tithi_radians( 30 , 1, 2026, 11.15, 5.5, &pdata);
//calculate_tithi_radians( 31 , 1, 2026, 8.416666666666666, 5.5, &pdata);
//calculate_tithi_radians( 1 , 2, 2026, 5.866666666666667, 5.5, &pdata);
//calculate_tithi_radians( 2 , 2, 2026, 3.6333333333333333, 5.5, &pdata);
//calculate_tithi_radians( 3 , 2, 2026, 1.8666666666666667, 5.5, &pdata);
//calculate_tithi_radians( 4 , 2, 2026, 0.6666666666666666, 5.5, &pdata);
//calculate_tithi_radians( 5 , 2, 2026, 0.15, 5.5, &pdata);
//calculate_tithi_radians( 6 , 2, 2026, 0.36666666666666664, 5.5, &pdata);
//calculate_tithi_radians( 7 , 2, 2026, 1.3, 5.5, &pdata);
//calculate_tithi_radians( 8 , 2, 2026, 2.9, 5.5, &pdata);
//calculate_tithi_radians( 9 , 2, 2026, 5.016666666666667, 5.5, &pdata);
//calculate_tithi_radians( 10 , 2, 2026, 7.45, 5.5, &pdata);
//calculate_tithi_radians( 11 , 2, 2026, 9.966666666666667, 5.5, &pdata);
//calculate_tithi_radians( 12 , 2, 2026, 12.366666666666667, 5.5, &pdata);
//calculate_tithi_radians( 13 , 2, 2026, 14.433333333333334, 5.5, &pdata);
//calculate_tithi_radians( 14 , 2, 2026, 16.016666666666666, 5.5, &pdata);
//calculate_tithi_radians( 15 , 2, 2026, 17.083333333333332, 5.5, &pdata);
//calculate_tithi_radians( 16 , 2, 2026, 17.566666666666666, 5.5, &pdata);
//calculate_tithi_radians( 17 , 2, 2026, 17.5, 5.5, &pdata);
//calculate_tithi_radians( 18 , 2, 2026, 16.95, 5.5, &pdata);
//calculate_tithi_radians( 19 , 2, 2026, 15.966666666666667, 5.5, &pdata);
//calculate_tithi_radians( 20 , 2, 2026, 14.633333333333333, 5.5, &pdata);
//calculate_tithi_radians( 21 , 2, 2026, 13.016666666666667, 5.5, &pdata);
//calculate_tithi_radians( 22 , 2, 2026, 11.166666666666666, 5.5, &pdata);
//calculate_tithi_radians( 23 , 2, 2026, 9.15, 5.5, &pdata);
//calculate_tithi_radians( 24 , 2, 2026, 7.033333333333333, 5.5, &pdata);
//calculate_tithi_radians( 25 , 2, 2026, 4.85, 5.5, &pdata);
//calculate_tithi_radians( 26 , 2, 2026, 2.6833333333333336, 5.5, &pdata);
//calculate_tithi_radians( 27 , 2, 2026, 0.55, 5.5, &pdata);
//calculate_tithi_radians( 27 , 2, 2026, 22.55, 5.5, &pdata);
//calculate_tithi_radians( 28 , 2, 2026, 20.716666666666665, 5.5, &pdata);
//calculate_tithi_radians( 1 , 3, 2026, 19.15, 5.5, &pdata);
//calculate_tithi_radians( 2 , 3, 2026, 17.916666666666668, 5.5, &pdata);
//calculate_tithi_radians( 3 , 3, 2026, 17.116666666666667, 5.5, &pdata);
//calculate_tithi_radians( 4 , 3, 2026, 16.816666666666666, 5.5, &pdata);
//calculate_tithi_radians( 5 , 3, 2026, 17.05, 5.5, &pdata);
//calculate_tithi_radians( 6 , 3, 2026, 17.883333333333333, 5.5, &pdata);
//calculate_tithi_radians( 7 , 3, 2026, 19.283333333333335, 5.5, &pdata);
//calculate_tithi_radians( 8 , 3, 2026, 21.183333333333334, 5.5, &pdata);
//calculate_tithi_radians( 9 , 3, 2026, 23.45, 5.5, &pdata);
//calculate_tithi_radians( 11 , 3, 2026, 1.9, 5.5, &pdata);
//calculate_tithi_radians( 12 , 3, 2026, 4.316666666666666, 5.5, &pdata);
//calculate_tithi_radians( 13 , 3, 2026, 6.483333333333333, 5.5, &pdata);
//calculate_tithi_radians( 14 , 3, 2026, 8.183333333333334, 5.5, &pdata);
//calculate_tithi_radians( 15 , 3, 2026, 9.266666666666667, 5.5, &pdata);
//calculate_tithi_radians( 16 , 3, 2026, 9.683333333333334, 5.5, &pdata);
//calculate_tithi_radians( 17 , 3, 2026, 9.383333333333333, 5.5, &pdata);
//calculate_tithi_radians( 18 , 3, 2026, 8.416666666666666, 5.5, &pdata);
//calculate_tithi_radians( 19 , 3, 2026, 6.883333333333333, 5.5, &pdata);
//calculate_tithi_radians( 20 , 3, 2026, 4.866666666666667, 5.5, &pdata);
//calculate_tithi_radians( 21 , 3, 2026, 2.5166666666666666, 5.5, &pdata);
//calculate_tithi_radians( 21 , 3, 2026, 23.933333333333334, 5.5, &pdata);
//calculate_tithi_radians( 22 , 3, 2026, 21.266666666666666, 5.5, &pdata);
//calculate_tithi_radians( 23 , 3, 2026, 18.633333333333333, 5.5, &pdata);
//calculate_tithi_radians( 24 , 3, 2026, 16.133333333333333, 5.5, &pdata);
//calculate_tithi_radians( 25 , 3, 2026, 13.833333333333334, 5.5, &pdata);
//calculate_tithi_radians( 26 , 3, 2026, 11.816666666666666, 5.5, &pdata);
//calculate_tithi_radians( 27 , 3, 2026, 10.116666666666667, 5.5, &pdata);
//calculate_tithi_radians( 28 , 3, 2026, 8.75, 5.5, &pdata);
//calculate_tithi_radians( 29 , 3, 2026, 7.766666666666667, 5.5, &pdata);
//calculate_tithi_radians( 30 , 3, 2026, 7.15, 5.5, &pdata);
//calculate_tithi_radians( 31 , 3, 2026, 6.933333333333334, 5.5, &pdata);
//calculate_tithi_radians( 1 , 4, 2026, 7.1, 5.5, &pdata);
//calculate_tithi_radians( 2 , 4, 2026, 7.683333333333334, 5.5, &pdata);
//calculate_tithi_radians( 3 , 4, 2026, 8.7, 5.5, &pdata);
//calculate_tithi_radians( 4 , 4, 2026, 10.15, 5.5, &pdata);
//calculate_tithi_radians( 5 , 4, 2026, 11.983333333333333, 5.5, &pdata);
//calculate_tithi_radians( 6 , 4, 2026, 14.166666666666666, 5.5, &pdata);
//calculate_tithi_radians( 7 , 4, 2026, 16.566666666666666, 5.5, &pdata);
//calculate_tithi_radians( 8 , 4, 2026, 19.016666666666666, 5.5, &pdata);
//calculate_tithi_radians( 9 , 4, 2026, 21.316666666666666, 5.5, &pdata);
//calculate_tithi_radians( 10 , 4, 2026, 23.25, 5.5, &pdata);
//calculate_tithi_radians( 12 , 4, 2026, 0.6166666666666667, 5.5, &pdata);
//calculate_tithi_radians( 13 , 4, 2026, 1.2833333333333332, 5.5, &pdata);
//calculate_tithi_radians( 14 , 4, 2026, 1.1333333333333333, 5.5, &pdata);
//calculate_tithi_radians( 15 , 4, 2026, 0.2, 5.5, &pdata);
//calculate_tithi_radians( 15 , 4, 2026, 22.516666666666666, 5.5, &pdata);
//calculate_tithi_radians( 16 , 4, 2026, 20.183333333333334, 5.5, &pdata);
}

