""" Use NexusAnalysis class and define \
     simplenexuscaller class which uses NexusAnalysis as a base data structure \
     and call peaks using that with a defined CLI same as chipr.
"""

import argparse
import sys
import numpy, pandas
from simplenexuscaller.nexusAnalysis import NexusAnalysis

class SimpleNexusCaller(object):
	""" Class for running the simplenexuscaller for peak calling on bedgraph \
		data from chip-nexus data.
	"""

	def __init__(self):
		""" Takes in command-line input from the user as chip-nexus bedgraph \
			 data.
		"""
		parser = argparse.ArgumentParser(prog='simplenexuscaller',
								description="Takes ChIP-nexus data in bedGraph "
											"format for + and - strand and performs "
											"fast and simple peak calling.\n")
		parser.add_argument("-i", "--input",
								help="ChIP-nexus bedGraph files, where each column is "
									 "[chr, start, end, count], with no column headers. "
									 "These files must be "
									 "in the order of counts on the + or - strand. "
									 "Inputs are separated by a single space. "
									 "Each BedGraph contains positions per base "
									 "where counts represent the number of 5' reads "
									 "mapping to that position. Where no reads "
									 "mapped, position refers to interval "
									 "where no reads mapped.",
								dest="input",
								type=str,
								nargs=2,
								required=True)
		parser.add_argument("-c", "--cutoff",
							help="Cutoff number of counts above which the position"
							"considered as having signal.",
							dest="cutoff",
							type=int,
							nargs=1,
							default=5,
							required=False)
		parser.add_argument("-f", "--falseInRowUpper",
							help="No. of no signal positions (count>cutoff) in row "
								 "before terminate extension of signal region.",
							dest="falseInRowUpper",
							type=int,
							nargs=1,
							default=10,
							required=False)
		parser.add_argument("-n", "--nInRowCutoff",
							help="The minimum length of a TF edge signal for the"
							"signal to be called as a true signal.",
							dest="nInRowCutoff",
							type=int,
							nargs=1,
							default=2,
							required=False)
		parser.add_argument("-d", "--distLimit",
							help="Minimum distance between signal range on the same"
							"strand for them to be considered the same or different"
							"TF binding signal edges.",
							dest="distLimit",
							type=int,
							nargs=1,
							default=40,
							required=False)
		parser.add_argument("-m", "--maxWidth",
							help="Maximum width a peak is allowed to be.",
							dest="maxWidth",
							type=int,
							nargs=1,
							default=100,
							required=False)
		parser.add_argument("-o", "--output",
							help="Output filename prefix. Automatically adds .bed",
							dest="output",
							type=str,
							default="simpleNexusPeaks",
							required=False)
		args = parser.parse_args(sys.argv[1:])
		self.runSimpleCaller(args)

	def runSimpleCaller(self, args):
		""" Defines how the simple caller runs based on user input.
		"""

		print("Reading in the data...")

		# Reading in the bedGraph data #
		posFileName, negFileName = args.input[0], args.input[1]
		colNames = ['chr', 'start', 'end', 'count']
		pos = pandas.read_csv(posFileName, sep='\t', names=colNames)
		neg = pandas.read_csv(negFileName, sep='\t', names=colNames)
		neg.loc[:, 'count'] = numpy.abs(neg.loc[:, 'count'])

		# Constructing the ChIP-nexus analysis object #
		nexus = NexusAnalysis(pos, neg)

		# Performing the peak calling #
		peaks = nexus.callPeaks(cutoff = args.cutoff,
								falseInRowUpper = args.falseInRowUpper,
								nInRowCutoff = args.nInRowCutoff,
								distLimit = args.distLimit,
								maxWidth = args.maxWidth)

		# Writing to file #
		nexus.write(f'{args.output}.bed')

def main():
	SimpleNexusCaller()

if __name__ == "__main__":
	main()








