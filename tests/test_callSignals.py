import unittest
from simplenexuscaller import callSignals
import numpy

class TestSignalFunctions(unittest.TestCase):

	def test_signalCalling(self):
		""" Tests whether the cuttoff based signal calling functioning as
		expected.
		"""
		counts = numpy.array( [1, 4, 5, 6, 8, 1, 5, 6] )
		cutoffs = [8, 5, 1]
		signals = callSignals.getSignalsInRange(counts, cutoffs, returnNumber=False)
		expected = [[False, False, False, False, True, False, False, False],
					[False, False, True, True, True, False, True, True],
					[True]*len(counts)]

		for i, signal in enumerate( signals ):
			expect = numpy.array( expected[i] )

			self.assertTrue( numpy.all(signal == expect) )

	def test_correctSignalTermination(self):
		""" Tests that when calling signalRanges terminates correctly
		"""
		signal1 = numpy.array( [False, False, True, True, True, False, False, False] )
		signal2 = numpy.array( [True, False, True, True, True, False, False, False] )
		signal3 = numpy.array( [True, False, True, True, True, False, False, True] )

		chrom1 = ['chr1']*len(signal1)
		chrom2 = ['chr1']+['chr2']* (len(signal1)-1)
		chrom3 = ['chr1']*3 + ['chr2']*5

		posLens = [1]*len(signal1)

		# Basic test
		ranges1 = callSignals.callSignalRanges(chrom1, posLens, signal1, 3, 2)
		self.assertTrue(ranges1 == [(2, 5)])

		# Test of minimum size constraint
		ranges1 = callSignals.callSignalRanges(chrom1, posLens, signal1, 4, 2)
		self.assertTrue(ranges1 == [])

		# Test false in row constraint
		ranges1 = callSignals.callSignalRanges(chrom1, posLens, signal2, 4, 2)
		self.assertTrue(ranges1 == [(0,5)])

		ranges1 = callSignals.callSignalRanges(chrom1, posLens, signal2, 4, 0)
		self.assertTrue(ranges1 == [])

		ranges1 = callSignals.callSignalRanges(chrom1, posLens, signal3, 4, 2)
		self.assertTrue(ranges1 == [(0, 8)])

		ranges1 = callSignals.callSignalRanges(chrom1, posLens, signal3, 4, 1)
		self.assertTrue(ranges1 == [(0, 5)])

		# Test the same chromosome constraint
		ranges1 = callSignals.callSignalRanges(chrom2, posLens, signal2, 3, 2)
		self.assertTrue(ranges1 == [(2, 5)])

		ranges1 = callSignals.callSignalRanges(chrom3, posLens, signal2, 3, 2)
		self.assertTrue(ranges1 == [(0, 3)])

		ranges1 = callSignals.callSignalRanges(chrom3, posLens, signal3, 3, 2)
		self.assertTrue(ranges1 == [(0, 3), (4, 8)])

		# Test for fact no counts represented with a range, whereas a count for
		# a position is represented with a count...
		posLens = [1, 100]+[1]*6 #case where one of the falses is very long.
		ranges1 = callSignals.callSignalRanges(chrom1, posLens, signal3, 4, 1)
		self.assertTrue(ranges1 == [])

		ranges1 = callSignals.callSignalRanges(chrom1, posLens, signal3, 3, 1)
		self.assertTrue(ranges1 == [(2, 5)])

		posLens = [1]*5 + [2, 2, 1] #case where two falses in a row add up
		ranges1 = callSignals.callSignalRanges(chrom1, posLens, signal3, 4, 4)
		self.assertTrue(ranges1 == [(0, 8)])

		ranges1 = callSignals.callSignalRanges(chrom1, posLens, signal3, 4, 3)
		self.assertTrue(ranges1 == [(0, 5)])

	def test_rangeSummits(self):
		""" Tests calling range summits.
		"""
		baseCounts = numpy.array( [1, 2, 3, 3, 3, 3, 5, 3, 10, 7] )
		ranges = [(0, 3), (3, 6), (4, 8), (7, 10)]
		vals = [baseCounts[ran[0]: ran[1]] for ran in ranges]
		summits = callSignals.getRangeSummits(baseCounts, ranges)

		for i, val in enumerate( vals ):
			self.assertTrue(val[summits[i]] == max(val))

if __name__ == '__main__':
	unittest.main()
