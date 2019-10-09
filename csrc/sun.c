#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>

#define floatingpoint double

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
	slon = sun(dd, mm, yy, zhr);
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

/*
 * MATLAB code from ESRL spreadsheet
---- https://ch.mathworks.com/matlabcentral/fileexchange/62180-sunriseset-lat-lng-utcoff-date-plot
function [sun_rise_set, varargout] = sunRiseSet( lat, lng, UTCoff, date, PLOT)
%SUNRISESET Compute apparent sunrise and sunset times in seconds.
%     sun_rise_set = sunRiseSet( lat, lng, UTCoff, date) Computes the *apparent* (refraction
%     corrected) sunrise  and sunset times in seconds from mignight and returns them as
%     sun_rise_set.  lat and lng are the latitude (+ to N) and longitude (+ to E), UTCoff is the
%     timezone, i.e. the local time offset to UTC (Coordinated Universal Time) in hours, and date is
%     the date in format 'dd-mmm-yyyy' ( see below for an example).
% 
%     [sun_rise_set, noon] = sunRiseSet( lat, lng, UTCoff, date) additionally returns the solar noon
%     in seconds from midnight.
% 
%     [sun_rise_set, noon, opt] = sunRiseSet( lat, lng, UTCoff, date) additionally returns the
%     information opt, which contains information on every second of the day:
%       opt.elev_ang_corr   : Apparent (refraction corrected) solar elevation in degrees
%       opt.azmt_ang        : Solar azimuthal angle (deg cw from N)
%       opt.solar_decl      : Solar declination in degrees
% 
%     sun_rise_set = sunRiseSet( ..., PLOT) If PLOT is true, plots of the elevation and azimuthal
%     angle are created.
% 
% EXAMPLE:
%     lat = 47.377037;    % Latitude (Zurich, CH)
%     lng = 8.553952;     % Longitude (Zurich, CH)
%     UTCoff = 2;         % UTC offset
%     date = '15-jun-2017';
% 
%     [sun_rise_set, noon, opt] = sunRiseSet( lat, lng, UTCoff, date, 1);
%
% 
% Richard Droste
% 
% Reverse engineered from the NOAA Excel:
% (https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html)
% 
% The formulas are from:
% Meeus, Jean H. Astronomical algorithms. Willmann-Bell, Incorporated, 1991.
% Process input
try
    nDays = daysact('30-dec-1899',  date);  % Number of days since beginning of the calculation
catch
    nDays = datenum(date)-datenum('30-dec-1899');
end
nTimes = 24*3600;                       % Number of seconds in the day
tArray = linspace(0,1,nTimes);
if nargin < 5
    PLOT = false;
end
% Compute
% Letters correspond to colums in the NOAA Excel
E = tArray;
F = nDays+2415018.5+E-UTCoff/24;
G = (F-2451545)/36525;
I = mod(280.46646+G.*(36000.76983+G*0.0003032),360);
J = 357.52911+G.*(35999.05029-0.0001537*G);
K = 0.016708634-G.*(0.000042037+0.0000001267*G);
L = sind(J).*(1.914602-G.*(0.004817+0.000014*G))+sind(2*J).* ...
    (0.019993-0.000101*G)+sind(3*J)*0.000289;
M = I+L;
P = M-0.00569-0.00478*sind(125.04-1934.136*G);
Q = 23+(26+((21.448-G.*(46.815+G.*(0.00059-G*0.001813))))/60)/60;
R = Q+0.00256*cosd(125.04-1934.136*G);
T = asind(sind(R).*sind(P));
U = tand(R/2).*tand(R/2);
V = 4*rad2deg(U.*sin(2*deg2rad(I))-2*K.*sin(deg2rad(J))+4*K.*U.*sin(deg2rad(J)).* ...
    cos(2*deg2rad(I))-0.5.*U.*U.*sin(4*deg2rad(I))-1.25.*K.*K.*sin(2.*deg2rad(J)));
AB = mod(E*1440+V+4*lng-60*UTCoff,1440);
if AB/4 < 0, AC = AB/4+180;else, AC = AB/4-180; end
AD = acosd(sind(lat)*sind(T)+cosd(lat)*cosd(T).*cosd(AC));
W = acosd(cosd(90.833)./(cosd(lat)*cosd(T))-tand(lat)*tand(T));
X = (720-4*lng-V+UTCoff*60)*60;
% Results in seconds
[~,noon]    = min(abs(X - nTimes*tArray));
[diff_sr, sunrise] = min(abs(X-round(W*4*60) - nTimes*tArray));
[diff_ss, sunset] = min(abs(X+round(W*4*60) - nTimes*tArray));
% Results in degrees
if nargout > 2 || PLOT
    solar_decl = T;
    elev_ang_corr = 90-AD;
    AC_ind = AC > 0;
    azmt_ang = mod(acosd(((sind(lat)*cosd(AD))-sind(T))./ ...
        (cosd(lat)*sind(AD)))+180,360);
    azmt_ang_2 = mod(540-acosd(((sind(lat)*cosd(AD))-sind(T))./ ...
        (cosd(lat)*sind(AD))),360);
    azmt_ang(~AC_ind) = azmt_ang_2(~AC_ind);
end
% Set sunrise and sunset values to -1 if out of 24h bounds
if abs(diff_sr) < 1
    sr_string = sprintf('Sunrise: %s', datestr(sunrise/nTimes,'HH:MM:SS'));
else
    sunrise = -1;
    sr_string = sprintf('Sunrise out of 24h bounds');
end
if abs(diff_ss) < 1
    ss_string = sprintf('Sunset: %s', datestr(sunset/nTimes,'HH:MM:SS'));
else
    sunset = -1;
    ss_string = sprintf('Sunset out of 24h bounds');
end
fprintf('%s\n%s\n', sr_string, ss_string)
% Generate output
sun_rise_set = [sunrise sunset];
if nargout > 1
    varargout{1} = noon;
end
if nargout > 2
    opt.elev_ang_corr = elev_ang_corr;
    opt.azmt_ang = azmt_ang;
    opt.solar_decl = solar_decl;
    varargout{2} = opt;
end
% Generate plots
if PLOT
    figure; hold on
    plot(linspace(0,24,nTimes), elev_ang_corr);
    xlabel('Hour'), ylabel('Angle [Deg]')
    xlim([0 24]), grid on
    title('Corrected Elevation Angle')
    
    figure;
    plot(linspace(0,24,nTimes), azmt_ang);
    xlabel('Hour'), ylabel('Angle [Deg]')
    xlim([0 24]), grid on
    title('Azimuthal Angle')
end
*/
