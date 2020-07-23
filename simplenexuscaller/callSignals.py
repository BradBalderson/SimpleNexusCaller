""" This script stores functions for calling ChIP-nexus signals.
Signals are defined as positions which could indicate an instance where
the edge of a TF bound to the DNA has been detected.

Signal ranges indicate
a range of such potential boundaries due to either technical variability
in 5' exonuclease activity, or biological variability of slightly different
TF binding position.

Variables for signal range calling include:

 1) Setting a cutoff for the number of counts at a position which is considered
 																	   a signal.
 2) Setting cutoff for how many signals in a row is minimum for a signal range.
 3) Setting cutoff for how many lack of signals in a row (by bp) results in
 	termination of a signal range.
"""

import numpy

def getSignalsInRange(counts, cuttoffs, returnNumber=True):
	""" Counts the number of bases which have a count greater than each of the \
	cutoffs. Has option to return either the number of bases or a boolean \
	numpy array which indicates where the base is greater than the cutoff.

	Args:
		counts (numpy.array<int>): Array indicating number of read counts \
									assigned to each position in the genome.

		cutoffs (list<int>): List of integers indicating different counts \
										above which is considered a signal.

		returnNumber (bool): Whether to return the number of signals detected \
									or as a boolean array indicating where \
									signals detected in the counts array.

	Returns:
		list<int> or list<numpy.array<bool>>: Former if returnNumber=True, \
											  latter if returnNumber=False. \
											Former is list indicating the \
											number of signals detected for each \
											cutoff. Latter indicates positions \
											where signals occur for each cutoff.
	"""

	nSignals = []
	for i, cutoff in enumerate( cuttoffs ):
		peakBool = counts>=cutoff

		if returnNumber:
			nSignals.append( len( numpy.where( peakBool )[0] ) )
		else:
			nSignals.append( peakBool )

	return nSignals

def getRangeSummits(counts, signalRanges):
	""" Calls the position with the largest signal in the signalRanges.

	Args:
		counts (numpy.array<int>): As in getSignalsInRange.

		signalRanges (list<tuple<int, int>>): As indicated in output of \
															callSignalRangesBed.

	Returns:
		list<int>: For each signal range, the number of positions from the \
					signal start position where the maximum count occurs; \
					such that signalRanges[0]+output[0] = position where max \
					count is in that signal range.
	"""

	signalSummits = []
	for (start, end) in signalRanges:
		maxSignal = max( counts[start:end] )
		signalSummit = numpy.where(counts[start:end]
								   == maxSignal)[0][0]
		signalSummits.append( signalSummit )

	return signalSummits

def callSignalRanges(posChroms, posLens, signals,
					 nInRowCutoff, falseInRowUpper):
	""" Calls signal ranges using method described in documentation for \
	callSignalRangesBed but using more basic input.

	Args:
		posChroms ( numpy.array<str> ): Specifies chromosome of each position.

		posLens ( numpy.array<int> ): Specifies length of each position.

		signals ( numpy.array<bool> ): Specifies which position considered to \
										have signal.

	Returns:
		list<tuple<int, int>>: signalRanges, as indicated in output of \
															callSignalRangesBed.
	"""

	signalLocs = numpy.where(signals)[0]

	# resolving multiple 'peaks' in a row by taking the one with the highest value
	signalRanges = [] #stores signal width
	for start in signalLocs:
		# Checking the start not already within a signal range
		if len(signalRanges) > 0 and \
				start >= signalRanges[-1][0] and start <= signalRanges[-1][1]:
			continue

		chrom = posChroms[start]
		end = start + 1
		falseInRow = 0
		falseInRowPosCounts = 0
		while falseInRow <= falseInRowUpper and end < len(signals) and \
			  chrom == posChroms[end]:

			if end not in signalLocs:
				falseInRow += posLens[end] #consider how long there is no signal
				falseInRowPosCounts += 1
			else:
				falseInRow = 0
				falseInRowPosCounts = 0

			end += 1


		nInRow = end - falseInRowPosCounts - start
		if nInRow >= nInRowCutoff:
			rangeEnd = end - falseInRowPosCounts
			signalRanges.append((start, rangeEnd))

	return signalRanges


def callSignalRangesBed(bedFrame, cutoff, falseInRowUpper, nInRowCutoff):
	"""Calls ChIP-nexus boundaries; i.e. the edge of where exonuclease stops \
	cutting. Imagine that Chip-Nexus peaks are really showing variation in the \
	boundary of the TF, either upstream (+) strand, or downstream boundary (-) \
	strand. This function aims to identify the range of boundaries on a \
	particular strand.

	Works by:

		1. Identifying bases which are greater than the cuttoff (signalLocs).

		2. Find the signalRanges: when find a signal, extend the \
		   signal downstream until get more than falseInRowUpper bases without \
		   signals. Also terminate if next bases refer to a different chromosome.

		3. Also find signalSummits within ranges: positions in the range which \
		   have the greatest number of counts. These are useful to find the \
		   most likely, or atleast the most frequent, TF boundary.

	Args:
		bedFrame (pandas.DataFrame) : Colnames are [chr, start, end, count]. \
									  Each row refers to a base with counts.

		cutoff (int): Upper cuttoff for a base to be considered as containing \
					  a signal.

		falseInRowUpper (int): Number of bases where no signal has occured in \
							   a row before terminating the signalRange.

		nInRowCutoff (int): Minimum number of signals in a row for to be \
							considered a signalRange.

	Returns:
		list<tuple<int, int>>, list<int>: List specifying the signalRanges, \
						which contains tuples referring to the startRow in the \
						inputted bedFrame where the signal begins, and the \
						endRow in the inputted bedFrame which is where the \
						signal ends. The second list contains the position \
						where the maximum signal occurs within the signal range \
						is zero based from the start position of the signal range.
						This is meant to indicate the most likely position for \
						the identified boundary.
	"""

	counts = bedFrame.values[:,3].astype(int) #Counts per position
	chroms = bedFrame.values[:,0].astype(str) #Chromosome of position
	#Lengths of each position
	posLens = bedFrame.values[:,2].astype(int) - bedFrame.values[:,1].astype(int)

	# Getting positions which have values greater than the cutoff
	signals = getSignalsInRange(counts, [cutoff], returnNumber=False)[0]

	# Calling the signalRanges
	signalRanges = callSignalRanges(chroms, posLens, signals,
									nInRowCutoff, falseInRowUpper)

	# Calling the summits
	signalSummits = getRangeSummits(counts, signalRanges)

	return signalRanges, signalSummits
