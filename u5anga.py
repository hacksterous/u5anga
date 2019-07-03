# (c) 2019 Anirban Banerjee
#Licensed under:
#GNU GENERAL PUBLIC LICENSE
#Version 3, 29 June 2007
from mpapbf import *
import gc
import mp5anga
from math import *
class u5anga ():
	pi = mpap('3.1415926535897932')
	pix2 = mpap('6.2831853071795865')
	D2R = pi/180
	R2D = mpap(180)/pi

	Longitude = mpap(81.889)
	Latitude = mpap(25.426)
	Elevation = mpap(0)

	Long = 81.889
	Lat = 25.426
	Elev = 0

	month = ["January","February","March","April","May","June",
	   "July","August","September","October","November","December"]

	rashi = ["Mesha","Vrisha","Mithuna","Karka","Simha","Kanya","Tula",
	   "Vrischika","Dhanu","Makara","Kumbha","Meena"]

	vaara = ["Ravi","Soma","Mangal","Budh","Brihaspati","Shukra","Shani"]

	tithi = ["Prathamaa","Dvitiya","Tritiya","Chaturthi","Panchami",
			"Shashthi","Saptami","Ashtami","Navami","Dashami","Ekadashi",
			"Dvadashi","Trayodashi","Chaturdashi","Purnima","Pratipada",
			"Dvitiya","Tritiya","Chaturthi","Panchami","Shashthi",
			"Saptami","Ashtami","Navami","Dashami","Ekadashi","Dvadashi",
			"Trayodashi","Chaturdashi","Amaavasya"]

	karana = ["Bava","Baalava","Kaulava","Taitula","Garija","Vanija",
	   "Vishti","Shakuni","Chatushpada","Naga","Kimstughna"]

	yoga = ["Vishakumbha","Priti","Ayushman","Saubhagya","Shobhana",
	   "Atiganda","Sukarman","Dhriti","Shula","Ganda","Vriddhi",
	   "Dhruva","Vyaghata","Harshana","Vajra","Siddhi","Vyatipata",
	   "Variyan","Parigha","Shiva","Siddha","Saadhya","Shubha","Shukla",
	   "Brahma","Indra","Vaidhriti"]

	nakshatra = ["Ashvini","Bharani","Krittika","Rohini","Mrigashira","Ardra",
			"Punarvasu","Pushya","Ashlesa","Magha","Purva Phalguni","Uttara Phalguni",
			"Hasta","Chitra","Svaati","Vishakha","Anuradha","Jyeshtha","Mula",
			"Purva Ashadha","Uttara Ashadha","Shravana","Dhanishtha","Shatabhisha",
			"Purva Bhaadra","Uttara Bhaadra","Revati"]

	def calc5 (self, dd, mm, yy, hr, zhr):
		mp5anga.set_date(dd, mm, yy)
		mp5anga.set_hour(hr, zhr)
		
		vaara = self.vaara[mp5anga.vaara()] 
		t = mp5anga.tithi()
		tithi = self.tithi[t]
		paksha = "Shukla" if (t <= 14) else "Krishna"
			
		nakshatra = self.nakshatra[mp5anga.nakshatra()]
		
		yoga = self.yoga[mp5anga.yoga()]
		
		karana = self.karana[mp5anga.karana()]
			
		#Calculate the rashi in which the moon is present
		rashi = self.rashi[mp5anga.rashi()]

		print(' Tithi: ' + tithi, '\n', paksha + ' paksha\n', vaara + 'vaara\n', 'Nakshatra: ' + nakshatra + '\n',\
				'Yoga: ' + yoga + '\n',\
				'Karana: ' + karana + '\n', 'Raashi: ' + rashi)
		return

	def sun (self, dd, mm, yy):
		#the Julian date
		t0 = (mpap(mm) - mpap(9)) // mpap(7)
		gc.collect()
		t1 = ((mpap(7)) * (mpap(yy) + mpap(5001) + t0)) // 4
		gc.collect()
		t2 = (mpap(275) * mpap(mm)) // mpap(9)
		d = mpap(367) * mpap(yy) - t1 
		gc.collect()
		d = d + t2 + mpap(dd) + mpap('1729777')
		#print ("Julian day to d is ", d)

		#d = Jdate is the Julian date
		#n is the number of days since Jan 1st, 2000 12:00.
		#2451545.0 is the equivalent Julian year of Julian days for 2000, 1, 1.5.
		#(68.184) / 86400 = 0.0008 is the fractional Julian Day for leap 
		#seconds and terrestrial time.
		n = d - mpap('2451545') + mpap('68.184') / mpap('86400')

		#---mean solar noon
		#jSTAR is an approximation of mean solar time at d,
		#expressed as a Julian day with the day fraction.
		#self.Longitude is the longitude west (west is negative, 
		#east is positive) of the observer on the Earth in degrees
		jSTAR = n - self.Longitude/mpap(360)
		#print ("Mean solar day is ", jSTAR)

		#---solar mean anomaly
		gc.collect()
		M = ((mpap('357.5291') + mpap('0.98560028') * jSTAR) % 360) * self.D2R

		#---equation of the centre
		#C is the Equation of the center value needed to calculate lambda (see next equation).
		#1.9148 is the coefficient of the Equation of the Center for the planet the observer is on (in this case, Earth)
		gc.collect()
		C = (mpap('1.9148') * M.sin() + mpap('0.02') * (M * mpap(2)).sin() + mpap('0.0003') * \
				(M * mpap(3)).sin()) * self.D2R
		#print ("equation of the centre is", C)

		#---ecliptic longitude -- lambda
		#lambdaLong is the ecliptic longitude.
		#102.9372 is a value for the argument of perihelion.
		gc.collect()
		lambdaLong = (M + C + self.pi + mpap('102.9372') * self.D2R) % self.pix2
		#print ("ecliptic longitude -- lambda is ", lambdaLong)

		#---solar transit
		#jTRANSIT is the Julian date for the local true solar transit (or solar noon).
		#2451545.0 is noon of the equivalent Julian year reference.
		#mpap(0.0053) * M.sin() - mpap('0.0069') * (lambdaLong * 2).sin() is a 
		#simplified version of the equation of time. The coefficients are fractional day minutes.
		gc.collect()
		jTRANSIT = mpap('2451545') + jSTAR + mpap(0.0053) * M.sin() - mpap('0.0069') * (lambdaLong * 2).sin()
		#print ("Julian date for the local true solar transit is ", jTRANSIT)

		#---declination of the sun
		#delta is the declination of the sun
		#23.44° is Earth's maximum axial tilt toward the sun
		gc.collect()
		sineDelta = lambdaLong.sin() * (mpap('23.44') * self.D2R).sin()
		gc.collect()
		delta = mpap(0) #sineDelta.asin()
		delta = sineDelta.asin()
		#print ("declination of the sun is", delta)

		#---elevation correction (elevation is in metres)
		#This corrects for both apparent dip and terrestrial refraction. 
		#For example, for an observer at 10,000 feet, add (−115°/60°) or about −1.92° to −0.83°.
		#gc.collect()
		elevationCorr = mpap('-2.076') / mpap(60)
		gc.collect()
		elevationCorr = elevationCorr * self.Elevation.sqrt() * self.D2R
		#---hour angle
		#omega is the hour angle from the observer's zenith;
		#phi is the north latitude of the observer (north is positive, 
		#south is negative) on the Earth.
		gc.collect()
		phi = self.Latitude * self.D2R
		#print("north latitude of the observer", phi)
		gc.collect()
		cosO = ((mpap('-0.83') * self.D2R + elevationCorr).sin() - phi.sin() * sineDelta)
		#print("cos0", cosO)
		gc.collect()
		cosOmega = cosO / (phi.cos() * delta.cos())
		#print ("cosine of hour angle from the observer's zenith is", cosOmega)
		gc.collect()
		omega = cosOmega.acos()
		#print ("hour angle from the observer's zenith is", omega)

		#calculate sunrise and sunset times
		gc.collect()
		jRISE = jTRANSIT - omega / self.pix2
		#print ("Rise Time: ", jRISE)
		gc.collect()
		jSET  = jTRANSIT + omega / self.pix2
		#print ("Set Time: ", jSET)

		gc.collect()
		jNOON = (jSET + jRISE) / mpap(2)
		gc.collect()
		dayLightHours = (jSET - jRISE) * mpap(24)
		#print ("dayLightHours is :", dayLightHours)
		gc.collect()
		riseToNoonHrs = (jNOON - jRISE) * mpap(24)
		gc.collect()
		noonToSetHrs = mpap(jSET - jNOON) * mpap(24)

		sriseHr = str( mpap(12) - riseToNoonHrs.ceil()) + ':' + ((mpap(1) - riseToNoonHrs.frac()) * 60).roundstr(2)
		ssetHr = str(int(noonToSetHrs) + 12) + ':' + (noonToSetHrs.frac() * 60).roundstr(2)
		print (" Sunrise at: ", sriseHr)
		print (" Sunset at: ", ssetHr)
		print (" Hours of daylight:", dayLightHours.roundstr(2))

		return

	def s2 (self, dd, mm, yy):
		#For debug
		d2r = 3.14159265/180
		r2d = 180/3.14159265
		#the Julian date
		t1 = (7 * (yy + 5001 + (mm - 9) // 7)) // 4
		t2 = (275 * mm) // 9
		d = 367 * yy - t1 + t2 + dd + 1729777
		#print ("Julian day to d is ", d)

		#d = Jdate is the Julian date
		#n is the number of days since Jan 1st, 2000 12:00.
		#2451545.0 is the equivalent Julian year of Julian days for 2000, 1, 1.5.
		#(68.184) / 86400 = 0.0008 is the fractional Julian Day for leap 
		#seconds and terrestrial time.
		n = d - 2451545.0 + (68.184/86400)
		#print ("n: ", n)

		#---mean solar noon
		#jSTAR is an approximation of mean solar time at d,
		#expressed as a Julian day with the day fraction.
		#self.Longitude is the longitude west (west is negative, 
		#east is positive) of the observer on the Earth in degrees
		jSTAR = n - self.Long/360.0
		#print ("jSTAR: ", jSTAR)

		#---solar mean anomaly
		M = ((357.5291 + 0.98560028 * jSTAR) % 360) * d2r
		#print ("M: ", M)

		#---equation of the centre
		#C is the Equation of the center value needed to calculate lambda (see next equation).
		#1.9148 is the coefficient of the Equation of the Center for the planet the observer is on (in this case, Earth)
		C = (1.9148 * sin(M) + 0.02 * sin(M*2) + 0.0003 * sin(M*3)) * d2r
		#print ("C: ", C)

		#---ecliptic longitude -- lambda
		#lambdaLong is the ecliptic longitude.
		#102.9372 is a value for the argument of perihelion.
		lambdaLong = (M + C + 3.14159265 + 102.9372 * d2r) % (2 * 3.14159265)
		#print ("lambdaLong: ", lambdaLong)

		#---solar transit
		#jTRANSIT is the Julian date for the local true solar transit (or solar noon).
		#2451545.0 is noon of the equivalent Julian year reference.
		#(0.0053) * sin(M) - (0.0069) * sim(lambdaLong * 2) is a 
		#simplified version of the equation of time. The coefficients are fractional day minutes.

		#jTRANSIT = 2451545.0 + jSTAR + 0.0053 * sin(M) - 0.0069 * sin(lambdaLong * 2)
		
		#make the noon of the equivalent Julian year reference as 0
		jTRANSIT = jSTAR + 0.0053 * sin(M) - 0.0069 * sin(lambdaLong * 2)
		
		#print ("jTRANSIT: ", jTRANSIT)

		#---declination of the sun
		#delta is the declination of the sun
		#23.44° is Earth's maximum axial tilt toward the sun
		sineDelta = sin(lambdaLong) * sin(23.44 * d2r)
		delta = asin(sineDelta)
		#print ("delta: ", delta)

		#---elevation correction (elevation is in metres)
		#This corrects for both apparent dip and terrestrial refraction. 
		#For example, for an observer at 10,000 feet, add (−115°/60°) or about −1.92° to −0.83°.
		elevationCorr = (-2.076 / 60) * sqrt(self.Elev) * d2r
		
		#---hour angle
		#omega is the hour angle from the observer's zenith;
		#phi is the north latitude of the observer (north is positive, 
		#south is negative) on the Earth.
		phi = self.Lat * d2r
		#print ("phi: ", phi)
		cosOmega = sin((-0.83 * d2r + elevationCorr) - sin(phi) * sineDelta) / \
					(cos(phi) * cos(delta))
		omega = acos(cosOmega)
		#print ("omega: ", omega)

		#calculate sunrise and sunset times
		jRISE = jTRANSIT - omega / (2 * 3.14159265)
		jSET  = jTRANSIT + omega / (2 * 3.14159265)

		#print ("Rise Time: ", jRISE)
		#print ("Set Time: ", jSET)

		jNOON = (jSET + jRISE) / 2
		dayLightHours = (jSET - jRISE) * 24
		#print ("jNOON is :", jNOON)
		riseToNoonHrs = (jNOON - jRISE) * 24
		noonToSetHrs = (jSET - jNOON) * 24

		sriseHr = str( 12 - ceil(riseToNoonHrs)) + ':' + str((1 - (riseToNoonHrs - int(riseToNoonHrs))) * 60)
		ssetHr = str(int(noonToSetHrs) + 12) + ':' + str((noonToSetHrs - int(noonToSetHrs)) * 60)
		print (" Sunrise at: ", sriseHr)
		print (" Sunset at: ", ssetHr)
		print (" Hours of daylight:", dayLightHours)

		return

#--- dd, mm, yy = 2,7,2019
#--- s2(2,7,2019)
#--- def s3 (dd, mm, yy):
#---     t1 = ((mpap(7)) * (mpap(yy) + mpap(5001) + (mpap(mm) - mpap(9)) // 7)) // 4
#---     t2 = (mpap(275) * mm) // 9
#---     d = mpap(367) * mpap(yy) - mpap(t1)
#--- 
#--- from u5anga import *; u=u5anga()
#--- u.sun(2,7,2019)
