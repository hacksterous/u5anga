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
		
	def tithis_to_file (self, startyr, filename, howmany):
		f = open(filename,"w")
		self.list_tithis(startyr, f, howmany)
		f.close()

	def all_tithis (self, d, m, y):
		#find all tithi changes starting at midnight on this date
		oldtithi = ""
		i = 0
		for minutes in range(1440):
			timescale = mp5anga.ts_at_mn (d, m, y) + (minutes/1440)
			#print ("timescale is ", timescale)
			tithi = self.tithis[mp5anga.tithi(timescale)]
			if tithi != oldtithi:
				newm = str(int(minutes%60))
				if len(newm) < 2:
					newm = "0"+newm
				if oldtithi != "":
					print (str(d)+"-"+str(m)+"-"+str(y)+" -> "+tithi+" starting "+str(int(minutes/60))+":"+newm)
				oldtithi = tithi
				i += 1

	def list_tithis (self, startyr, howmany=1, fobj=None):
		daylist = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
		i = 0
		oldtithi = ""
		while True:
			for m in range(12):
				for d in range(daylist[m]):
					leapyr = (startyr % 4 == 0) and ((startyr % 100 != 0) or (startyr % 400 == 0))
					for minutes in range(1440):
						timescale = mp5anga.ts_at_mn (d+1, m+1, startyr) + (minutes/1440)
						tithi = self.tithis[mp5anga.tithi(timescale)]
						if tithi != oldtithi:
							if oldtithi != "AA":
								newm = str(int(minutes%60))
								if len(newm) < 2:
									newm = "0"+newm
								print (str(d+1)+"-"+str(m+1)+"-"+str(startyr)+","+str(int(minutes/60))+":"+newm+","+tithi)
								if fobj != None:
									fobj.write (str(d+1)+"-"+str(m+1)+"-"+str(startyr)+","+str(int(minutes/60))+":"+newm+","+tithi+"\n")
							i += 1
							if i >= howmany:
								return
							oldtithi = tithi
					if m == 1 and d == 27 and not leapyr:
						break
					elif m == 1 and d == 28 and leapyr:
						break
					elif m == 11 and d == 30:
						startyr += 1
						break

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
		if sun[2] == "0":
			print ("No sunrise today.")
		elif sun[2] == "24":
			print ("No sunset today.")
		else:
			print (" Sunrise at: ", sun[0])
			print (" Sunset at : ", sun[1])
			print (" Daylight hours:", sun[2])


