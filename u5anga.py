# (c) 2019 Anirban Banerjee
#from u5anga import *; u = u5anga(5.5, 12.93340, 77.59630)
#Licensed under:
#GNU GENERAL PUBLIC LICENSE
#Version 3, 29 June 2007
from mpapbf import *
import gc
import mp5anga
from math import *
gc.enable()

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

	vaaras = ["Ravi","Soma","Mangal","Budh","Brihaspati","Shukra","Shani"]

	tithis = ["Pratipad","Dvitiya","Tritiya","Chaturthi","Panchami",
			"Shashthi","Saptami","Ashtami","Navami","Dashami","Ekadashi",
			"Dvadashi","Trayodashi","Chaturdashi","Purnima","Pratipad",
			"Dvitiya","Tritiya","Chaturthi","Panchami","Shashthi",
			"Saptami","Ashtami","Navami","Dashami","Ekadashi","Dvadashi",
			"Trayodashi","Chaturdashi","Amaavasya"]

	karana = ["Bava","Baalava","Kaulava","Taitula","Garija","Vanija",
	   "Vishti","Shakuni","Chatushpada","Naga","Kimstughna"]

	yogas = ["Vishakumbha","Preeti","Ayushman","Saubhagya","Shobhana",
	   "Atiganda","Sukarman","Dhriti","Shula","Ganda","Vriddhi",
	   "Dhruva","Vyaghata","Harshana","Vajra","Siddhi","Vyatipata",
	   "Variyan","Parigha","Shiva","Siddha","Saadhya","Shubha","Shukla",
	   "Brahma","Indra","Vaidhriti"]

	nakshatras = ["Ashvini","Bharani","Krittika","Rohini","Mrigashira","Ardra",
			"Punarvasu","Pushya","Ashlesa","Magha","Purva Phalguni","Uttara Phalguni",
			"Hasta","Chitra","Svaati","Vishakha","Anuradha","Jyeshtha","Mula",
			"Purva Ashadha","Uttara Ashadha","Shravana","Dhanishtha","Shatabhisha",
			"Purva Bhaadra","Uttara Bhaadra","Revati"]

	def __init__ (self, zhr, latt, longt):
		mp5anga.set_zone(zhr, latt, longt)

	def mrashi (self, dd, mm, yyyy, hr):
		print (self.rashi[mp5anga.vaara(mp5anga.ts_at_mn (dd, mm, yyyy) + (hr/24.0))])

	def vaara (self, dd, mm, yyyy):
		print (self.vaaras[mp5anga.vaara(mp5anga.ts_at_mn (dd, mm, yyyy))])

	def ayanansha (self, dd, mm, yyyy, hr):
		print ("ayanansha =", mp5anga.ayanansha(mp5anga.ts_at_mn (dd, mm, yyyy) + (hr/24.0)))

	def nakshatra (self, dd, mm, yyyy, hr):
		print ("nakshatra =", self.nakshatra[mp5anga.nakshatra(mp5anga.ts_at_mn (dd, mm, yyyy) + (hr/24.0))])
		
	def tithi (self, dd, mm, yyyy, hr):
		timescale = mp5anga.ts_at_mn (dd, mm, yyyy) + (hr/24.0)
		t = mp5anga.tithi(timescale)
		#print ("tithi index t =", t)
		print(' Tithi: ' + ("Shukla " if (t <= 14) else "Krishna ") + self.tithis[t])

	def sun (self, dd, mm, yyyy):
		#print ("sriseHr: ", sriseHr)
		timescale = mp5anga.ts_at_mn (dd, mm, yyyy) + 2451544 #2451543.5 + 12/24 (at noon)
		sun = mp5anga.sun(timescale).split(' ')
		#print ("mpap sriseHr: ", sriseHr)
		print (" Sunrise at: ", sun[0])
		print (" Sunset at : ", sun[1])
		print (" Daylight hours:", sun[2])


