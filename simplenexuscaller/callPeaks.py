""" Based on tf binding site boundaries identified on the + and - strand,
identifies likely binding sites by linking + strand boundaries with their closest
downstream negative strand boundary. Any resulting binding sites (peaks) with
a distance greater than a particular cutoff are excluded, since these are
called as being far too long to be reasonably considered a binding site for
the TF.
"""

import numpy, pandas

#TODO sanity check where mix these two around... should give same results.
def getCandidatePeaks(posBoundaries, negBoundaries):
	""" For each boundary on the + strand, match to closest downstream boundary \
	on - strand. Note that the positions of the boundaries are specified in a \
	particular order, this allows for the first negative boundary downstream \
	of the positive boundary on the same chromosome can be taken as the \
	closest matching boundary.

	Args:
		posBoundaries (pandas.DataFrame): Specify locations of a called tf \
										binding boundary on the + strand.
										Colnames are \
										[chr, start, end, count, originIndex]. \
							Column 'originIndex' which specifies the row index \
					where that position can be found in the original bedFrame.

		negBoundaries (pandas.DataFrame): Same as posBoundaries, except \
									  specifies TF binding boundary on - strand.

	Returns:
		list<list<str, int, int, int, int, int>>): List of called tf binding \
						events 'peaks'. Each peak specified by: \
						 [chrom, start, end, width, originIndex1, originIndex2].
	"""
	peaks = []
	lastNegi = 0 #Keeps track of the last position detected downstream of +
				 #strand position.
	for posi in range(posBoundaries.shape[0]):
		posBoundary = posBoundaries.values[posi, :] #+ strand tf boundary
		# posStart = posBoundary[1]
		#
		# # Getting positions downstream...
		# negStarts = negBoundaries.loc[:,'start'].values
		# negIndicesDownstream = numpy.where(negStarts > posStart)[0]

		for negi in range(lastNegi, negBoundaries.shape[0]): #negIndicesDownstream:
			negBoundary = negBoundaries.values[negi, :] #- strand tf boundary

			# if on the same chromosome and - boundary downstream of +
			chromEqual = posBoundary[0] == negBoundary[0]
			isDownstream = negBoundary[1] > posBoundary[1]
			if chromEqual and isDownstream:
				# chrom, start, end, width, originIndex1, originIndex2
				peak = [posBoundary[0], posBoundary[1], negBoundary[1],
						negBoundary[1]-posBoundary[1],
						posBoundary[-1], negBoundary[-1]]
				peaks.append(peak)
				lastNegi = negi
				break

	return numpy.array(peaks)

def getPeaks(posBoundaries, negBoundaries, maxWidth):
	""" Calls peaks by matching tf binding boundaries on + strand with closest \
	boundary on - strand, and filtering these based on a minimum width.

	Args:
		posBoundaries (pandas.DataFrame): As indicated in 'getCandidatePeaks'.
		negBoundaries (pandas.DataFrame): As indicated in 'getCandidatePeaks'.
		maxWidth (int): Maximum width a peak is allowed to be.

	Returns:
		pandas.DataFrame: Dataframe of called tf binding events \
							('peaks'). Each peak specified by per row as: \
							[chrom, start, end, originIndex1, originIndex2].
	"""

	# Getting candidate peaks #
	peaks = getCandidatePeaks(posBoundaries, negBoundaries)

	# Getting candidates which meet minWidth threshold #
	widths = peaks[:, 3].transpose().astype(int)
	widthBool = widths < maxWidth
	peakBed = pandas.DataFrame(peaks[widthBool,:],
							   columns=['chr', 'start', 'end', 'width',
										'originIndex1', 'originIndex2'])

	return peakBed
