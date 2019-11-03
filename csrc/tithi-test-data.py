import sys
with open("tithi-test-data2.csv") as fp:
	line = fp.readline()
	while line:
		f = line.strip().split(',')
		try:
			if "AM" in f[5] and int(f[3]) == 12:
				f[3] = "0"
			if "PM" in f[5] and int(f[3]) < 12:
				f[3] = str(int(f[3]) + 12 + int(f[4])/60)
			else:
				f[3] = str(int(f[3]) + int(f[4])/60)
			print ("calculate_tithi_radians(" + f[1] + ", " + str(int(f[0])) + ", " + str(int(f[2]) + 2000) + ", " + f[3]  + ", 5.5, &pdata);")
			line = fp.readline()
		except:
			print ("Error with line ", line.strip())
			sys.exit()
