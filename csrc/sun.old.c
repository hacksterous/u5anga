//#########################################################
//sunrise function using LibBF math library
//#########################################################


#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>

/*
#include "libbf.h"

#ifdef BARE_M
#include "mpconfig.h"
#include "misc.h"
#define realloc m_realloc
#define free m_free
#endif

static bf_t A, B, R, S;
static bf_context_t bf_ctx;

static void *my_bf_realloc(void *opaque, void *ptr, size_t size) {
    return realloc(ptr, size);
}

void bf_initialize (void) {
    bf_context_init(&bf_ctx, my_bf_realloc, NULL);
}

void bf_exit (void) {
    bf_context_end(&bf_ctx);
}
*/

#define floatingpoint float

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

floatingpoint sun (floatingpoint* sunriseLSTp, 
			floatingpoint* sunsetLSTp, 
			int dd, 
			int mm, 
			int yyyy, 
			floatingpoint zhr, 
			floatingpoint latt, 
			floatingpoint longt) {

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

/*
char* sunbf (const char* dd, const char* mm, const char* yyyy) {

	const uint32_t precisiondigits = 15;
	const uint32_t rnd_mode = 0;
    int32_t prec;
    uint32_t status;
    char *digits;
	bf_t t1, t2, d, n, temp1, temp2, temp3;
	char *D2RSTR = '17.4532925199432957692369';

    prec = (limb_t) (precisiondigits * 4 + 32);

    bf_init(&bf_ctx, &t1);
    bf_init(&bf_ctx, &t2);
    bf_init(&bf_ctx, &d);
    bf_init(&bf_ctx, &n);
    bf_init(&bf_ctx, &temp1);
    bf_init(&bf_ctx, &temp2);
    bf_init(&bf_ctx, &temp3);



    bf_ftoa(&digits, &result, 10, precisiondigits + 1,
                             BF_FTOA_FORMAT_FIXED | rnd_mode);

    bf_delete(&t);
    bf_delete(&o);
    bf_delete(&l);
    bf_delete(&temp1);
    bf_delete(&temp2);
    bf_delete(&temp3);
    bf_delete(&ayan);
    bf_delete(&D2R);

	return digits;
}

*/

int main () {
	floatingpoint srise, sset;
	floatingpoint sunhrs = sun (&srise, &sset, 
					6, 10, 2019, 
					//5.5, 65.972, 77.595);
					5.5, 66, 77.595);

	printf ("srise = %f\n", srise);
	printf ("sset = %f\n", sset);
	printf ("sunhrs = %f\n", sunhrs);

}

/*
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
