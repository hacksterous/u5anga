# (c) 2019 Anirban Banerjee
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

	def calc5 (self, dd, mm, yyyy, hr, zhr):
		mp5anga.set_date(dd, mm, yyyy)
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

	def sun (self, dd, mm, yyyy, zhr, latt, longt):
		if abs(latt) > 65:
			print ("Error: cannot do latitudes beyond the Arctic/Antarctic Circles")
			return 0.0
		mp5anga.set_date(dd, mm, yyyy)
		mp5anga.set_hour(-1, zhr) #set only time zone
		#print ("sriseHr: ", sriseHr)
		sriseHr = mpap(mp5anga.srise (latt, longt))
		#print ("mpap sriseHr: ", sriseHr)
		ssetHr = mpap(mp5anga.sset (latt, longt))
		#print ("sun sets: ", ssetHr)
		dayLightHours = mpap(mp5anga.shrs (latt, longt) / 60)
		#print ("sun rises fracx60: ", sriseHr.frac()*60)
		minutes = (sriseHr.frac()*60).roundstr(0)
		#print ("minutes: ", minutes)
		if len(minutes) == 1:
			minutes = '0' + minutes
		#print ("sun rises fracx60 roundstr0: ", (sriseHr.frac()*60).roundstr(0))
		print (" Sunrise at: ", (str(sriseHr.floor()+1) if minutes == "60" else str(sriseHr.floor())) + ":" + ("00" if minutes == "60" else minutes))
		minutes = (ssetHr.frac()*60).roundstr(0)
		if len(minutes) == 1:
			minutes = '0' + minutes
		print (" Sunset at : ", (str(ssetHr.floor()+1) if minutes == "60" else str(ssetHr.floor())) + ":" + ("00" if minutes == "60" else minutes))
		minutes = (dayLightHours.frac()*60).roundstr(0)
		if len(minutes) == 1:
			minutes = '0' + minutes
		print (" Daylight hours:", str(dayLightHours.floor())+":"+ ("00" if minutes == "60" else minutes) +" minutes")

		return sriseHr.float()

