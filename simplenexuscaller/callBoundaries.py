""" This script stores functions which use signal ranges and signal summits \
(derived from callSignals.py) and calls a boundary on where the most likely \
(or atleast most frequent) position where the edge of the TF binding occurs \
for each signal range.

It also resolves 'dual' boundaries, which are cases where two signal ranges \
occur very close to one another, indicating they are not likely different \
boundaries or are quite possibly two different boundaries.
"""

import numpy, pandas

def getBoundaries(signalRanges, signalSummits, bedFrame):
	""" Gets the positions where the signal range signal is largest, \
	which is indicative of a potential TF binding event boundary.

	Args:
		signalRanges ( list<tuple<int, int>> ): List specifying the
						signalRanges, contains tuples referring to the \
						startRow in the inputted bedFrame where the signal \
						begins, and the endRow in the inputted bedFrame which \
						is where the signal ends.

		signalSummits ( list<int> ): Contains the position where the maximum \
						signal occurs within the signal range. Is zero based \
						from the start position of the signal range. \
						This is meant to indicate the most likely position for \
						the identified boundary.

		bedFrame ( pandas.DataFrame ): Colnames are [chr, start, end, count]. \
									  Each row refers to a base with counts.

	Returns:
		pandas.DataFrame: bedFrame but subsetted to just where the summits of \
						the inputted signalRanges occur. Also includes an extra \
						column 'originIndex' which specifies the row index \
						where that position can be found in the original bedFrame.
	"""

	boundaries = []
	for i, (start, end) in enumerate(signalRanges):
		boundaries.append(start + signalSummits[i])

	boundaryFrame = bedFrame.iloc[boundaries, :]
	boundaryFrame.loc[:,'originIndex'] = list(boundaryFrame.index)

	return boundaryFrame

# TODO impliment the other options.
# TODO consider may represent different binding events from cellular \
#  heterogeneity in the sample which may result in a slightly different binding \
#  position.
def resolveDualBoundaries(boundaries, strand, distLimit, method):
	""" Resolves the occurence of dual boundaries which are neighbouring to one \
		another. Three methods are available:

		1) 'largestSignal' takes the boundary with the largest number of counts.

		2) 'wide' takes the boundary which is 5' most for + strand, and 3' \
					most for - strand. (Currently not implimented)

		3) 'narrow' takes the boundary which is 3' most for + strand, and 5' \
					most for - strand. (Currently not implimented)

	Args:
		boundaries ( pandas.DataFrame ): Bed format, rows are positions \
										specifying a boundary as a single base. \
										Column names are chromsome, start, end, \
										count, originIndex.

		strand ( str ): Whether boundary refers to + or - strand.

		distLimit ( int ): The distance between boundaries for them to be \
						   considered 'dual'.

		method ( str ): Either 'largestSignal', 'wide', 'narrow'; description \
						above.

	Returns:
		pandas.DataFrame: Same as the boundaries, except any dual boundaries \
						  removed so now only one boundary remains.
	"""

	if method == 'largestSignal':
		return resolveBoundariesWithSignal(boundaries, distLimit)


def resolveBoundariesWithSignal(boundaries, distLimit):
	""" Resolves the occurence of dual boundaries which are neighbouring to \
		one another by taking the peak with the highest number of counts.

	Args:
		boundaries ( pandas.DataFrame ): Bed format, rows are positions \
										specifying a boundary as a single base. \
										Column names are chromsome, start, end \
										count.

		distLimit ( int ): The distance between boundaries for them to be \
						   considered 'dual'.

	Returns:
		pandas.DataFrame: Same as the boundaries, except any dual boundaries \
						  removed so now only one boundary remains.
	"""

	# Resolving dual boundaries by taking the one with largest counts
	boundaryValues = boundaries.values
	collapsedBoundaries = []
	rowi = 0
	while rowi < boundaries.shape[0]-1:
		pos1 = boundaryValues[rowi, :]
		pos2 = boundaryValues[rowi + 1, :]

		# If two boundaries are within a certain range from one another
		chromEqual = pos1[0] == pos2[0] #chromosomes are same
		withinRange = pos2[1] - pos1[1] <= distLimit
		if chromEqual and withinRange:
			posCounts = numpy.array([pos1[2], pos2[2]])
			maxPos = numpy.where(posCounts == max(posCounts))[0][0]
			collapsedBoundaries.append( [pos1, pos2][maxPos] ) #add boundary
			rowi += 2 # skip over the next position

		else:
			collapsedBoundaries.append( pos1 ) #adding in the boundary
			rowi += 1 # go to the next position

	collapsedBoundaries = pandas.DataFrame(collapsedBoundaries,
									  				 columns=boundaries.columns)

	return collapsedBoundaries
