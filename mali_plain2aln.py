import sys, os, string

if __name__ == "__main__":
	# Length of each alignment segment
	LINE_WIDTH = 60

	# Get file data from stdin
	aliLines = sys.stdin.read().split("\n")
	while "" in aliLines: aliLines.remove("")
	while len(aliLines) > 2: del aliLines[2]

	# Sort data and get alignment length
	qArray = aliLines[0].split("\t")
	queryStart, queryAli, queryEnd, queryId = int(qArray[0]), qArray[1], int(qArray[2]), qArray[3]
	queryAli = queryAli.replace("-",".")
	sArray = aliLines[1].split("\t")
	sbjctStart, sbjctAli, sbjctEnd, sbjctId = int(sArray[0]), sArray[1], int(sArray[2]), sArray[3]
	sbjctAli = sbjctAli.replace("-",".")
	length = len(queryAli)

	CURR_POS = 0
	qFrom, qTo, sFrom, sTo = queryStart, 0, sbjctStart, 0
	while (CURR_POS < length):
		if length - CURR_POS < 60: LINE_WIDTH = length - CURR_POS
		qGaps, sGaps = 0, 0
		# Get alignment slices
		qAli = queryAli[CURR_POS:CURR_POS+LINE_WIDTH]
		sAli = sbjctAli[CURR_POS:CURR_POS+LINE_WIDTH]
		# Count gaps
		for x in range(0,len(qAli),1):
			if qAli[x] == ".": qGaps += 1
			if sAli[x] == ".": sGaps += 1
		qTo = qFrom + (LINE_WIDTH-1) - qGaps
		sTo = sFrom + (LINE_WIDTH-1) - sGaps
		print "\n"
		print "%-20s %-10d%s\t%d" %(queryId, qFrom, qAli, qTo)
		print "%-20s %-10d%s\t%d" %(sbjctId, sFrom, sAli, sTo)

		qFrom = qTo+1
		sFrom = sTo+1
		CURR_POS += LINE_WIDTH

	print "\n"
	sys.exit(0)
