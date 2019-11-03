import sys
import re
with open("tithi-input.txt") as fp:
	line = fp.readline().strip()
	while line:
		#try:
		if True:
			print (line)
			#r = re.compile('^\d+\/\d+\/\d+ \d+:\d+ [A|P]M\s+\d+\.\S+\s*')
			r = re.compile('^\d+\/\d+\/\d+ \d+:\d.*')
			line = fp.readline().strip()
			m = r.match(line)
			print (m)
		#except:
		else:
			print ("Error with line ", line)
			sys.exit()
		
