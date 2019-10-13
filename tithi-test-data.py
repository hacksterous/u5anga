with open("tithi-test-data.csv") as fp:
	line = fp.readline()
	while line:
		f = line.strip().split(',')
		if "AM" in f[5] and int(f[3]) == 12:
			f[3] = "0"
		if "PM" in f[5] and int(f[3]) < 12:
			f[3] = str(int(f[3]) + 12 + int(f[4])/60)
		else:
			f[3] = str(int(f[3]) + int(f[4])/60)
		print ("calculate_tithi (" + f[1] + ", " + str(int(f[0])) + ", " + str(int(f[2]) + 2000) + ", " + f[3]  + ", 5.5, &pdata);")
		line = fp.readline()
