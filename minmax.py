import sys
import getopt
import os.path

outputFilePath = None
writeMode = 'a'
freqdict = dict() #dictionary mapping codons to their frequencies
mapdict = dict() #dictionary mapping codons to amino acid
aatoFreq = dict() #dictionary mapping each amino acid to a list of the frequencies of possible codons for that amino acid
minMaxValues = [] #list holding the min/max value result for each sliding window
windowSize = 17 #size of the sliding window
helpMessage = '\nUsage: python minmax.py <frequencies file> <sequence file> [options]\n' \
			  'Use option -h for help.'
longHelp = '\nUsage: python minmax.py <frequencies file> <sequence file> [options]\n' \
		   '\nCalculates the minmax values for a codon sequence. Outputs a list with the minmax value for each window, with the first and last floor(windowSize/2) elements populated by "None".\n\n' \
		   '   frequencies file : A text file containing the individual codon appearance frequencies and associated amino acids for a particular species.\n' \
		   '   sequence file    : A text file containing the base sequence to be evaluated.\n\n' \
		   'Options:\n' \
		   '   -h      : Display the help screen (this).\n' \
		   '   -o      : Data will be written to "output.txt" instead of printed to the screen.\n' \
		   '   -O file : Data will be written to file instead of printed to the screen. If file does not exist, it will be created. Output will be appended to file\'s existing content.\n' \
		   '   -r      : When writing output to a file, replace (overwrite) the existing content instead of appending.\n' \
		   '   -w size : Sets the size of the sliding window (measured in codons). Must be between 2 and the total number of codons in the sequence. Defaults to 17.\n'

#parse command line arguments
try:
	opts, args = getopt.gnu_getopt(sys.argv[1:], "hoO:rw:")
except getopt.GetoptError:
	print(helpMessage)
	sys.exit(2)

#parse options
for opt, val in opts:
	if opt == '-h':
		print(longHelp)
		sys.exit()
	elif opt == '-o':
		outputFilePath = 'output.txt'
	elif opt == '-O':
		outputFilePath = val
	elif opt == '-r':
		writeMode = 'w'
	elif opt == '-w':
		try:
			num = int(val)
			if num < 2:
				print('Error: The window size must be at least 2.')
				sys.exit(2)
			else:
				windowSize = num
		except ValueError:
			print('Error: The window size must be an integer.')
			sys.exit(2)

#check to see if the required input files are provided
if len(args) < 2:
	print(helpMessage)
	sys.exit(2)

#resolve the input files
if os.path.isfile(args[0]):
	frequenciesFilePath = args[0]
else:
	print('Error: Frequencies file not found.')
	sys.exit(2)

if os.path.isfile(args[1]):
	sequenceFilePath = args[1]
else:
	print('Error: Sequence file not found.')
	sys.exit(2)

frequenciesFile = list(open(frequenciesFilePath))

#parse frequencies files for relevant information
for line in frequenciesFile:
	parts = line.split()
	freqdict[parts[0]] = float(parts[2])
	mapdict[parts[0]] = parts[1]
	aatoFreq[parts[1]] = []

#Populate the first floor(windowsize/2) entries of minMaxValues with None
for i in range(int(windowSize/2)):
	minMaxValues.append(None)

#Adds frequencies for each amino acid to its corresponding list in the dictionary
for line in frequenciesFile:
	parts = line.split()
	aatoFreq[parts[1]].append(float(parts[2]))

#For a given input fasta sequence (will chose an arbitrary one)
sequenceFile = open(sequenceFilePath)
codons = []
codstr = ""
for line in sequenceFile:
	line = line.rstrip()
	codstr += line
i = 0
while (i+3) <= len(codstr):
		codons.append(codstr[i:i+3])
		i += 3

if windowSize > len(codons):
	print('Error: the input sequence contains', len(codons), 'codons.')
	print('The window size must not be larger than the number of codons in the sequence.')
	sys.exit(2)

#Using the specified sliding window size (floor(windowSize/2)) on either side of the central codon), min/max is calculated
for i in range(len(codons)-windowSize+1):
	window = codons[i:i+windowSize] #list of the codons in the current window

	Actual = 0.0     #average of the actual codon frequencies
	Max = 0.0        #average of the min codon frequencies
	Min = 0.0        #average of the max codon frequencies
	Avg = 0.0        #average of the averages of all the frequencies associated with each amino acid

	#Sum the frequencies
	for codon in window:
		frequencies = aatoFreq[mapdict[codon]] #list of all frequencies associated with the amino acid this codon encodes

		Actual += freqdict[codon]
		Max += max(frequencies)
		Min += min(frequencies)
		Avg += sum(frequencies)/len(frequencies)

	#Divide by the window size to get the averages
	Actual = Actual/windowSize
	Max = Max/windowSize
	Min = Min/windowSize
	Avg = Avg/windowSize

	percentMax = ((Actual-Avg)/(Max-Avg))*100
	percentMin = ((Avg-Actual)/(Avg-Min))*100

	if percentMax >= 0:
		minMaxValues.append(percentMax)
	else:
		minMaxValues.append(-percentMin)

#Populate the last floor(windowsize/2) entries of minMaxValues with None
for i in range(int(windowSize/2)):
	minMaxValues.append(None)

outputText = 'Frequencies file: ' + frequenciesFilePath + '\n' + \
			 'Sequence file: ' + sequenceFilePath + '\n' + \
			 'Number of codons in sequence: ' + str(len(codons)) + '\n' + \
	    	 'Sliding window size: ' + str(windowSize) + '\n' + \
			 'Number of windows evaluated: ' + str(len(codons) - windowSize + 1)
print('\n' + outputText)

#print the values or write them to a file, as requested
if outputFilePath is None:
	print('\nMin Max Values:')
	print(minMaxValues)
else:
	print('Data written to "' + outputFilePath + '"\n')
	outputFile = open(outputFilePath, writeMode)
	outputFile.write(outputText + '\n\nMin Max Values:\n')
	for value in minMaxValues:
		outputFile.write(str(value))
		outputFile.write('\n')
	outputFile.write('\n\n')
	outputFile.close()
