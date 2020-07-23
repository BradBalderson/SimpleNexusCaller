from simplenexuscaller import callSignals, callBoundaries, callPeaks

class NexusAnalysis(object):
	""" A datastructure for holding nexus bedGraph data and performing \
	the analysis on them by using callSignals, callBoundaries, and callPeaks.
	Stores all of the intermediate states that the data passes through on way \
	to peak calling.

	Construction is by contract, no error checking.
	"""

	pos = None
	neg = None

	def __init__(self, pos, neg):
		""" NexusAnalysis object constructor.

			Args:
				pos (pandas.DataFrame):
								Bedgraph of ChIP-nexus counts on the positive \
								strand. Ascending order of the genome \
								positions. Columns are ['chr', 'start', 'end', \
								'count']. The 'counts' column refers to the \
								number of reads which had their 5' end start \
								at that particular base. Where the 'count' \
								column is non-zero for a particular row, will \
								refer to one specific base. Where count is \
								zero for a particular row, will refer to a large \
								region of bases where there are non-zero counts \
								until the next position where a count occurs.

				neg (pandas.DataFrame): Same as pos, except on the negative strand.
		"""
		self.pos = pos
		self.neg = neg

	def callPeaks(self, cutoff=10, falseInRowUpper=10, nInRowCutoff=2,
				  distLimit=40, maxWidth=100):
		""" Performs peak calling on ChIP-nexus data.

		Args:
			cutoff (int): No. of counts to be considered signal

			falseInRowUpper (int): No. of no signals in row before terminate.

			nInRowCutoff (int): Min no. of signals in row to be considered \
																    signalRange.

			distLimit (int): The distance between boundaries for them to be \
						   considered 'dual'; i.e. possible that there is \
						   two binding sites right next to one another, or \
						   could be due to techical noise.

			maxWidth (int): Maximum width a peak is allowed to be. \
							Recommend to set this about 2-3 times the TF \
							binding site width.
		"""

		print("Calling TF signals...")
		# Calling 'signals' (defined as positions which could indicate an
		# instance where the edge of a TF bound to the DNA has been detected.)
		self.signalRanges, self.signalSummits = \
			callSignals.callSignalRangesBed(self.pos,
										  cutoff, falseInRowUpper, nInRowCutoff)

		self.signalRangesNeg, self.signalSummitsNeg = \
			callSignals.callSignalRangesBed(self.neg,
										  cutoff, falseInRowUpper, nInRowCutoff)

		print("Calling TF binding boundaries...")
		# Calling the 'boundaries' (where the most likely \
		# (or atleast most frequent) position where the edge of the
		# TF binding occurs for each signal range.)
		self.posBounds = callBoundaries.getBoundaries(self.signalRanges,
												 	  self.signalSummits,
												 	  self.pos)
		self.negBounds = callBoundaries.getBoundaries(self.signalRangesNeg,
												 	  self.signalSummitsNeg,
												 	  self.neg)

		print("Resolving dual boundaries...\n")
		self.posBoundaries = callBoundaries.resolveBoundariesWithSignal(
																self.posBounds,
																distLimit)
		self.negBoundaries = callBoundaries.resolveBoundariesWithSignal(
																self.negBounds,
																distLimit)
		print("TF boundaries detected on + and - strand (respectively):")
		print(self.posBoundaries.shape[0], self.negBoundaries.shape[0])
		print("Numbers should be roughly the same if chosen parameters are good"
			  " for the dataset.\n")

		print("Calling peaks...\n")
		# Call peaks by matching tf binding boundaries on + strand with closest
		# boundary on - strand, and filtering these based on a minimum width.
		self.peaks = callPeaks.getPeaks(self.posBoundaries, self.negBoundaries,
										maxWidth)
		print(f"Detected {len(self.peaks)} peaks.")

		return self.peaks

	def write(self, fileName):
		""" Writes the peaks to a bed file with columns: chr, start, end.
		"""
		if type(self.peaks) == type(None):
			print("Need to call callPeaks() first.")
			return

		peakBedFile = self.peaks.loc[:, ['chr', 'start', 'end']]
		peakBedFile.to_csv(fileName, sep='\t', index=False, header=False)


