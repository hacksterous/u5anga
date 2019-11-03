import math
def timescale(dd, mm, yy, hr, zhr):
	#https:#stjarnhimlen.se/comp/ppcomp.html#3
	#The time scale in these formulae are counted in days. 
	#Hours, minutes, seconds are expressed as fractions of a day. 
	#Day 0.0 occurs at 2000 Jan 0.0 UT (or 1999 Dec 31, 24:00 UT). 
	#This "day number" d is computed as follows (y=year, m=month, D=date, 
	#UT=UT in hours+decimals):
	d = 367*yy \
			- 7 * ( yy + (mm+9)/12 ) / 4 \
			- 3 * ( ( yy + (mm-9)/7 ) / 100 + 1 ) / 4 \
			+ 275*mm/9 + dd - 730515
	return (d + (hr - zhr)/24.0)#add 2451543.5 to get Julian Date;


def radians (degrees):
	return (3.14159265358979323846 * degrees / 180.0)


def degrees (radians):
	return (180.0 * radians /3.14159265358979323846)

def sun_esrl (julian_day, zhr, latt, longt):

	#Reverse engineered from the NOAA Excel:
	#https:#www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html

	julian_century = (julian_day - 2451545.0)/36525.0 #G
	#printf ("julian_century is %.15f\n", julian_century)

	sun_geom_mean_long_deg = (280.46646+julian_century* \
							(36000.76983 + julian_century*0.0003032) % 360.0) #I
	#printf ("sun_geom_mean_long_deg is %.15f\n", sun_geom_mean_long_deg)

	sun_geom_mean_anom_deg = 357.52911+julian_century* \
							(35999.05029 - 0.0001537*julian_century) #J
	#printf ("sun_geom_mean_anom_deg is %.15f\n", sun_geom_mean_anom_deg)

	earth_orbit_eccentricity = 0.016708634-julian_century* \
							(0.000042037+0.0000001267*julian_century) #K
	#printf ("earth_orbit_eccentricity is %.15f\n", earth_orbit_eccentricity)

	sun_eqn_of_centre = math.sin(radians(sun_geom_mean_anom_deg))* \
							(1.914602-julian_century*(0.004817+0.000014*julian_century))+ \
							math.sin(radians(2*sun_geom_mean_anom_deg))* \
							(0.019993-0.000101*julian_century)+ \
							math.sin(radians(3*sun_geom_mean_anom_deg))*0.000289 #L
	#printf ("sun_eqn_of_centre is %.15f\n", sun_eqn_of_centre)

	sun_true_long_deg = sun_geom_mean_long_deg + sun_eqn_of_centre #M
	#printf ("sun_true_long_deg is %.15f\n", rev(sun_true_long_deg))

	sun_true_anom_deg = sun_geom_mean_anom_deg + sun_eqn_of_centre #N
	#printf ("sun_true_anom_deg is %.15f\n", sun_true_anom_deg)

	sun_rad_vector_au = (1.000001018*(1-earth_orbit_eccentricity*earth_orbit_eccentricity))/ \
							(1+earth_orbit_eccentricity*math.cos(radians(sun_true_anom_deg))) #O
	#printf ("sun_rad_vector_au is %.15f\n", sun_rad_vector_au)

	sun_apparent_long_deg = sun_true_long_deg-0.00569-0.00478* \
							math.sin(radians(125.04-1934.136*julian_century)) #P
	#printf ("sun_apparent_long_deg is %.15f\n", rev(sun_apparent_long_deg))

	mean_oblique_ecliptic_deg = 23+(26+((21.448-julian_century*(46.815+ \
							julian_century*(0.00059-julian_century*0.001813))))/60.0)/60.0 #Q
	#printf ("mean_oblique_ecliptic_deg is %.15f\n", mean_oblique_ecliptic_deg)

	oblique_correction_deg = mean_oblique_ecliptic_deg+ \
							0.00256*math.cos(radians(125.04-1934.136*julian_century)) #R
	#printf ("oblique_correction_deg is %.15f\n", oblique_correction_deg)

	#printf ("radians(sun_apparent_long_deg) is %.15f\n", radians(sun_apparent_long_deg))

	sun_right_ascension_deg = degrees(math.atan2(math.cos(radians(oblique_correction_deg))*math.sin(radians(sun_apparent_long_deg)),
							math.cos(radians(sun_apparent_long_deg)))) #S
	#printf ("sun_right_ascension_deg is %.15f\n", sun_right_ascension_deg)

	sun_declination_deg = degrees(math.asin(math.sin(radians(oblique_correction_deg))* \
							math.sin(radians(sun_apparent_long_deg)))) #T
	#print ("sun_declination_deg is %.15f\n" % sun_declination_deg)

	var_y = math.tan(radians(oblique_correction_deg/2))*math.tan(radians(oblique_correction_deg/2)) #U

	equation_of_time_min = 4*degrees(var_y*math.sin(2*radians(sun_geom_mean_long_deg))- \
							2*earth_orbit_eccentricity*math.sin(radians(sun_geom_mean_anom_deg))+ \
							4*earth_orbit_eccentricity*var_y*math.sin(radians(sun_geom_mean_anom_deg))* \
							math.cos(2*radians(sun_geom_mean_long_deg))- \
							0.5*var_y*var_y*math.sin(4*radians(sun_geom_mean_long_deg))- \
							1.25*earth_orbit_eccentricity*earth_orbit_eccentricity* \
							math.sin(2*radians(sun_geom_mean_anom_deg))) #V
	#printf ("equation_of_time_min is %.15f\n", equation_of_time_min)
	
	hour_angle_cosine = math.cos(radians(90.833))/(math.cos(radians(latt))* \
							math.cos(radians(sun_declination_deg)))-math.tan(radians(latt))* \
							math.tan(radians(sun_declination_deg))
	print ("hour_angle_cosine is ", hour_angle_cosine)

	if hour_angle_cosine > 1.0:
		#sun never rises
		return (0, 0, 0)
	elif (hour_angle_cosine < -1.0):
		#sun never sets
		return (0, 0, 24)

	hour_angle_sunrise_deg = degrees(math.acos(math.cos(radians(90.833))/(math.cos(radians(latt))* \
							math.cos(radians(sun_declination_deg)))-math.tan(radians(latt))* \
							math.tan(radians(sun_declination_deg)))) #W
	#printf ("hour_angle_sunrise_deg is %.15f\n", hour_angle_sunrise_deg)

	solar_noon_LST = (720.0-4*longt-equation_of_time_min+zhr*60)/1440.0 #X
	#printf ("solar_noon_LST is %.15f\n", solar_noon_LST)

	sunriseLST = (solar_noon_LST - (hour_angle_sunrise_deg*4)/1440.0) * 24; #Y
	sunsetLST = (solar_noon_LST + (hour_angle_sunrise_deg*4)/1440.0)*24; #Z
	sunlight_duration_hrs = hour_angle_sunrise_deg * 8 / 60 #AA

	return (sunriseLST, sunsetLST, sunlight_duration_hrs)


julian_day =  timescale (15, 10, 2019, 12.0, 5.5) + 2451543.5;
print (sun_esrl (julian_day, 5.5, 82.93340, 77.59630))
