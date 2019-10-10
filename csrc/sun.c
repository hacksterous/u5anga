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

struct panchanga {
	char dtoday[16];
	char dvaara[16];
	char dyoga[64];
	char dnakshatra[64];
	char dtithi[64];
	char dkarana[64];
	char dpaksha[64];
	char drashi[32];
};

#define D2R (M_PI/180.0)
#define R2D (180.0/M_PI)
#define REVOLUTIONS(x)	((x)-floor((x)/360.0)*360.0)

#define floatingpoint double

static char month[][15]={"January","February","March","April","May","June",
		  "July","August","September","October","November","December"};

static char rashi[][15]={"Mesha","Vrishabha","Mithuna","Karka","Simha","Kanya","Tula",
		   "Vrischika","Dhanu","Makara","Kumbha","Meena"};

static char day[][15]={"Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday"};

static char tithi[][15]={"Prathamaa","Dvitiya","Trithiya","Chaturthi","Panchami",
		  "Shashthi","Saptami","Ashtami","Navami","Dashami","Ekadashi",
		  "Dvadashi","Trayodashi","Chaturdashi","Purnima","Pratipada",
		  "Dvitiya","Tritiya","Chaturthi","Panchami","Shashthi",
		  "Saptami","Ashtami","Navami","Dashami","Ekadashi","Dvadashi",
		  "Trayodashi","Chaturdashi","Amaavasya"};

static char karan[][15]={"Bava","Baalava","Kaulava","Taitula","Garija","Vanija",
		  "Vishti","Shakuni","Chatushpada","Naga","Kimstughna"};

static char yoga[][15]={"Vishakumbha","Priti","Ayushman","Saubhagya","Shobhana",
		 "Atiganda","Sukarman","Dhriti","Shula","Ganda","Vriddhi",
		 "Dhruva","Vyaghata","Harshana","Vajra","Siddhi","Vyatipata",
		 "Variyan","Parigha","Shiva","Siddha","Saadhya","Shubha","Shukla",
		 "Brahma","Indra","Vaidhriti"};

static char nakshatra[][20]={"Ashvini","Bharani","Krittika","Rohini","Mrigashira","Ardra",
		      "Punarvasu","Pushya","Ashlesa","Magha","Purva Phalguni","Uttara Phalguni",
		      "Hasta","Chitra","Svaati","Vishakha","Anuradha","Jyeshtha","Mula",
		      "Purva Ashadha","Uttara Ashadha","Shravana","Dhanishtha","Shatabhisha",
		      "Purva Bhaadra","Uttara Bhaadra","Revati"};


floatingpoint datenum(	int day, 
				int mon, 
				int year, 
				int hour, 
				int min, 
				int sec) {
	int tmp1, tmp2, tmp3;
	floatingpoint	tmp4, tmp5;
	floatingpoint dNum;
	int cumudays[] = {0, 0,31,59,90,120,151,181,212,243,273,304,334};
	/* Calculate the serial date number:*/
	tmp1 = 365 * year  + cumudays[mon] + day;
	tmp2 = year / 4 - year / 100 + year / 400;
	tmp3 = (year % 4 != 0) - (year % 100 != 0) + (year % 400 != 0);
	tmp4 = (floatingpoint) (tmp1+tmp2+tmp3);
	tmp5 = (hour * 3600 + min * 60 + sec) / 86400.0;
	dNum = tmp4 + tmp5;
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

double JulianDay (int date, int month, int year, double hr, double zhr) {
	//http://www.geoastro.de/elevazmoon/basics/meeus.htm
	double UT = hr - zhr;
    if (month<=2) {month=month+12; year=year-1;}
    return (int)(365.25*year) + (int)(30.6001*(month+1)) - 15.0 + 1720996.5 + date + UT/24.0;
}

floatingpoint JulianDay2 (int dd, 
			int mm, 
			int yyyy, 
			double hr,
			floatingpoint zhr) {

	unsigned int day = datenum (dd, mm, yyyy, 12, 0, 0) - 
							datenum (30, 12, 1899, 12, 0, 0);
	printf ("JulianDay2 is %d\n", day);

	const floatingpoint E = 0.5; //always calculate at noon = 12 hours from local midnight

	floatingpoint julian_day = day + 2415018.5 + E - (zhr/24.0); //F
	printf ("JulianDay2 is %f\n", julian_day);

}

floatingpoint sun_esrl (floatingpoint* sunriseLSTp, 
			floatingpoint* sunsetLSTp, 
			int dd, 
			int mm, 
			int yyyy, 
			floatingpoint zhr, 
			floatingpoint latt, 
			floatingpoint longt) {

	// Reverse engineered from the NOAA Excel:
	//(https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html)

	unsigned int day = datenum (dd, mm, yyyy, 12, 0, 0) - 
							datenum (30, 12, 1899, 12, 0, 0);
	printf ("day is %d\n", day);

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
	printf ("sun_true_long_deg is %f\n", sun_true_long_deg);

	floatingpoint sun_true_anom_deg = sun_geom_mean_anom_deg + sun_eqn_of_centre; //N
	//printf ("sun_true_anom_deg is %f\n", sun_true_anom_deg);

	floatingpoint sun_rad_vector_au = (1.000001018*(1-earth_orbit_eccentricity*earth_orbit_eccentricity))/
							(1+earth_orbit_eccentricity*cos(radians(sun_true_anom_deg))); //O
	//printf ("sun_rad_vector_au is %f\n", sun_rad_vector_au);

	floatingpoint sun_apparent_long_deg = sun_true_long_deg-0.00569-0.00478*
							sin(radians(125.04-1934.136*julian_century)); //P
	printf ("sun_apparent_long_deg is %f\n", sun_apparent_long_deg);

	floatingpoint mean_oblique_ecliptic_deg = 23+(26+((21.448-julian_century*(46.815+
							julian_century*(0.00059-julian_century*0.001813))))/60.0)/60.0; //Q
	//printf ("mean_oblique_ecliptic_deg is %f\n", mean_oblique_ecliptic_deg);

	floatingpoint oblique_correction_deg = mean_oblique_ecliptic_deg+
							0.00256*cos(radians(125.04-1934.136*julian_century)); //R
	//printf ("oblique_correction_deg is %f\n", oblique_correction_deg);

	//printf ("radians(sun_apparent_long_deg) is %f\n", radians(sun_apparent_long_deg));

	floatingpoint sun_right_ascension_deg = degrees(atan2(cos(radians(oblique_correction_deg))*sin(radians(sun_apparent_long_deg)),
							cos(radians(sun_apparent_long_deg)))); //S
	//printf ("sun_right_ascension_deg is %f\n", sun_right_ascension_deg);

	floatingpoint sun_declination_deg = degrees(asin(sin(radians(oblique_correction_deg))*
							sin(radians(sun_apparent_long_deg)))); //T
	printf ("sun_declination_deg is %f\n", sun_declination_deg);

	floatingpoint var_y = tan(radians(oblique_correction_deg/2))*tan(radians(oblique_correction_deg/2)); //U

	floatingpoint equation_of_time_min = 4*degrees(var_y*sin(2*radians(sun_geom_mean_long_deg))-
							2*earth_orbit_eccentricity*sin(radians(sun_geom_mean_anom_deg))+
							4*earth_orbit_eccentricity*var_y*sin(radians(sun_geom_mean_anom_deg))*
							cos(2*radians(sun_geom_mean_long_deg))-
							0.5*var_y*var_y*sin(4*radians(sun_geom_mean_long_deg))-
							1.25*earth_orbit_eccentricity*earth_orbit_eccentricity*
							sin(2*radians(sun_geom_mean_anom_deg))); //V
	//printf ("equation_of_time_min is %f\n", equation_of_time_min);
	
	float A = cos(radians(90.833));
	float B = cos(radians(latt))*cos(radians(sun_declination_deg));
	float C = tan(radians(latt))*tan(radians(sun_declination_deg));

	printf ("A = cos(radians(90.833) = %f\n", cos(radians(90.833)));
	printf ("B = cos(radians(latt))*cos(radians(sun_declination_deg)) = %f\n", cos(radians(latt))*cos(radians(sun_declination_deg)));
	printf ("A/B = %f\n", A/B);
	printf ("tan(radians(latt) = %f\n", tan(radians(latt)));
	printf ("tan(radians(sun_declination_deg) = %f\n", tan(radians(sun_declination_deg)));
	printf ("C = tan(radians(latt))*tan(radians(sun_declination_deg)) = %f\n", tan(radians(latt))*tan(radians(sun_declination_deg)));
	printf ("A/B - C = %f\n", A/B - C);
	floatingpoint hour_angle_sunrise_deg = degrees(acos(cos(radians(90.833))/(cos(radians(latt))*
							cos(radians(sun_declination_deg)))-tan(radians(latt))*
							tan(radians(sun_declination_deg)))); //W
	printf ("hour_angle_sunrise_deg is %f\n", hour_angle_sunrise_deg);

	floatingpoint solar_noon_LST = (720.0-4*longt-equation_of_time_min+zhr*60)/1440.0; //X
	//printf ("solar_noon_LST is %f\n", solar_noon_LST);

	*sunriseLSTp = (solar_noon_LST - (hour_angle_sunrise_deg*4)/1440.0) * 24; //Y

	*sunsetLSTp = (solar_noon_LST + (hour_angle_sunrise_deg*4)/1440.0)*24; //Z

	floatingpoint sunlight_duration_min = hour_angle_sunrise_deg * 8; //AA

	return (sunlight_duration_min);
}

double moon_long2(int dd, int mm, int yy, double hr, double zhr) {
	//https://astronomy.stackexchange.com/questions/24859/local-sidereal-time
	//https://astronomy.stackexchange.com/questions/30355/
	// how-to-calculate-the-ground-track-of-the-moons-position-on-the-earths-surface
	double right_ascension_moon; //equals local mean sidereal time (LMST)
	double d = JulianDay(dd, mm, yy, hr, zhr); //julian day
	double GMST = 100.4606184 + 0.9856473662862 * d + 15.0 * (hr - zhr);//greenwich MST
	return 0.0;
}

double fpart (double x) {
	double fx = x - floor(x);
	if (fx < 0)
		fx += 1;
	return fx;
}

static double moon_long(int dd, int mm, int yy, double hr, double zhr) {
	double l;
	return l;
}

static double sun_long(int dd, int mm, int yy, double hr, double zhr) {
	double l;
	return l;
}

double moon_ra (double t) {
//http://www.stargazing.net/kepler/riset.bas
// returns ra (and optionally dec) of Moon to 5 arc min (ra) and 1 arc min (dec)
// for a few centuries either side of J2000.0
// Predicts rise and set times to within minutes for about 500 years
// in past - TDT and UT time diference may become significant for long
// times
	double p2 = 6.283185307;
	double ARC = 206264.8062;
	double COSEPS = 0.91748;
	double SINEPS = 0.39778;
	double L0 = fpart(.606433 + 1336.855225 * t);     //mean long Moon in revs
	double L = p2 * fpart(.374897 + 1325.55241 * t);  //mean anomaly of Moon
	double LS = p2 * fpart(.993133 + 99.997361 * t);  //mean anomaly of Sun
	double d = p2 * fpart(.827361 + 1236.853086 * t); //diff longitude sun and moon
	double F = p2 * fpart(.259086 + 1342.227825 * t); //mean arg latitude
	//longitude correction terms
	double dL = 22640 * sin(L) - 4586 * sin(L - 2 * d);
	dL = dL + 2370 * sin(2 * d) + 769 * sin(2 * L);
	dL = dL - 668 * sin(LS) - 412 * sin(2 * F);
	dL = dL - 212 * sin(2 * L - 2 * d) - 206 * sin(L + LS - 2 * d);
	dL = dL + 192 * sin(L + 2 * d) - 165 * sin(LS - 2 * d);
	dL = dL - 125 * sin(d) - 110 * sin(L + LS);
	dL = dL + 148 * sin(L - LS) - 55 * sin(2 * F - 2 * d);
	//latitude arguments
	double S = F + (dL + 412 * sin(2 * F) + 541 * sin(LS)) / ARC;
	double h = F - 2 * d;
	//latitude correction terms
	double N = -526 * sin(h) + 44 * sin(L + h) - 31 * sin(h - L) - 23 * sin(LS + h);
	N = N + 11 * sin(h - LS) - 25 * sin(F - 2 * L) + 21 * sin(F - L);
	double lmoon = p2 * fpart(L0 + dL / 1296000.0);  //Lat in rads
	double bmoon = (18520.0 * sin(S) + N) / ARC;     //long in rads
	//convert to equatorial coords using a fixed ecliptic
	double CB = cos(bmoon);
	double x = CB * cos(lmoon);
	double V = CB * sin(lmoon);
	double W = sin(bmoon);
	double y = COSEPS * V - SINEPS * W;
	double Z = SINEPS * V + COSEPS * W;
	double rho = sqrt(1.0 - Z * Z);
	double dec = (360.0 / p2) * atan(Z / rho);
	double ra = (48.0 / p2) * atan(y / (x + rho));
	if (ra < 0)
	    ra = ra + 24.0;
	return ra;
}

void calculate_tithi(int dd, int mm, int yy, double hr, double zhr, struct panchanga *pdata) {
	double slon, mlon, tmlon, tslon;
	int n;

	double d = JulianDay(dd, mm, yy, hr, zhr) -  JulianDay(31, 12, 1999, 0, 0);
	printf ("----\n");
	printf ("New Julian Day is %f\n", d);
	printf ("----\n");

	//Calculate, moon and sun longitude
	slon = sun_long(dd, mm, yy, zhr);
	printf ("sun long is %f\n", slon);
	mlon = moon_long(d);
	printf ("moon long is %f\n", mlon);
	mlon = mlon+((mlon<slon)?360:0);
	printf ("updated moon long is %f\n", mlon);
	//Calculate Tithi and Paksha
	tmlon = mlon+((mlon<slon)?360:0);
	tslon = slon;
	n = (int)((tmlon-tslon)/12);
	strcpy(pdata->dtithi,tithi[n]);
	(n<=14)?strcpy(pdata->dpaksha,"Shukla"):strcpy(pdata->dpaksha,"Krishna");
}

int main () {
	floatingpoint srise, sset;
	floatingpoint sunhrs = sun_esrl (&srise, &sset, 
					6, 10, 2019, 
					//5.5, 65.972, 77.595);
					5.5, 66, 77.595);

	printf ("srise = %f\n", srise);
	printf ("sset = %f\n", sset);
	printf ("sunhrs = %f\n", sunhrs);



}

