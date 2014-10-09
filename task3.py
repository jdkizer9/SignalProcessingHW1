import scipy
import scipy.io.wavfile
import pylab
import matplotlib
import numpy as np
import itertools
from scipy.signal import butter, lfilter 

kWindowSize = 20
kWindowLength = 4096
kWindowShift = 2048
kSamplingRate = 44100
kHighCut = 20000
kLowCut = 200

# Computes the Short-Time Fourier Transform (STFT) of a signal, with a given
# window length, and shift between adjacent windows
def stft(x, window_len=kWindowLength, window_shift=kWindowShift):
	w = scipy.hamming(window_len)
	X = scipy.array([scipy.fft(w*x[i:i+window_len])
		for i in range(0, len(x)-window_len, window_shift)])

	return np.array(scipy.absolute(X[:,0:window_len/2]))

# Plot a transformed signal X, i.e., if X = stft(x), then
# plot_transform(X) plots the spectrogram of x
def plot_transform(X):
	pylab.ion()
	pylab.figure()
	pylab.imshow(scipy.log(X.T), origin='lower', aspect='auto', interpolation='nearest', norm=matplotlib.colors.Normalize())
	pylab.xlabel('Window index')
	pylab.ylabel('Transform coefficient')
	pylab.ioff()


#t and f are indices
def isMaxima(t, f, X):
	Xdims = X.shape

	windowMinT = max(0, t-(kWindowSize/2))
	windowMinF = max(0, f-(kWindowSize/2))

	windowMaxT = min(Xdims[0], t+(kWindowSize/2))
	windowMaxF = min(Xdims[1], f+(kWindowSize/2))

	window = X[windowMinT:windowMaxT, windowMinF:windowMaxF]
	#slice array (inclusive)

	maxIndices = np.unravel_index(np.argmax(window),window.shape)
	if( (t == (maxIndices[0]+windowMinT)) and (f == (maxIndices[1]+windowMinF))):
		return True
	else:
		return False

def setOfPeaksForTimeFrequencyDataWithMagnitude(X):
	setOfMaxima = set()
	Xdims = X.shape
	for i in range(0, Xdims[0]):
		for j in range(0, Xdims[1]):
			if isMaxima(i, j, X):
				magnitude = X[i][j]
				setOfMaxima.add((i,j, magnitude))
	return setOfMaxima


# Plot a list of peaks in the form [(s1, f1), (s2, f2), ...]
def plot_peaks(peak_list,start,end):
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(1,1,1)    

    filteredPeakList = list(filter(lambda x: (x[0]>=start and x[0]<=end), peak_list))

    #s_list, f_list = zip(*filteredPeakList)

    #s_list, f_list = zip(*peak_list)
    s_list, f_list, m_list = zip(*peak_list)
    #print (s_list)



    matplotlib.pyplot.plot(s_list, f_list, 'bo',)    
    ymin, ymax = ax.get_ylim()    
    ax.vlines((start,end),ymin,ymax,'red')
    matplotlib.pyplot.xlabel('Window index')
    matplotlib.pyplot.ylabel('Transform coefficient')


def generate2DArrayOfFrequencyMagnitudePairs(mySetOfPoints, numWindows):
	array = [[] for i in range(numWindows)]
	for s in mySetOfPoints:
		array[s[0]].append((s[1], s[2]))
	return array

 
#map functions

#return the largest magnitude in (f, m) tuple array x
def largestMagnitude(x):
	returnMagnitude = 0.0
	#print(x)
	for t in x:
		#print(t)
		if t[1] > returnMagnitude: 
			returnMagnitude = t[1]
	return returnMagnitude

#return the largest magnitude in (f, m) tuple array x 
#w is the weighting function
#def largestWeightedMagnitude(w, x):


#return the highest frequency in (f, m) tuple array x
def highestFrequencyIndex(x):
	returnFrequency = 0.0
	for t in x:
		if t[0] > returnFrequency: 
			returnFrequency = t[0]
	return returnFrequency

#return the lowest frequency in (f, m) tuple array x
def lowestFrequencyIndex(x):
	if len(x) == 0:
		return 0.0
	else:
		returnFrequency = x[0][0]
		for t in x:
			if t[0] < returnFrequency: 
				returnFrequency = t[0]
		return returnFrequency

#return the sum of magnitudes in (f, m) tuple array x
def sumOfMagnitudes(x):
	magnitudes = [t[1] for t in x]
	return sum(magnitudes)

#return the weighted sum of magnitudes in (f, m) tuple array x
#def weightedSumOfMagnitues(w, x):

def crossCorrelate(x1, x2, numberOfWindows, mappingFunction):
	#print(x1)
	xArray1 = generate2DArrayOfFrequencyMagnitudePairs(x1, numberOfWindows)
	y1 = list(map(mappingFunction, xArray1))
	#print(len(y1))
	#print(y1)
	#print(x2)
	xArray2 = generate2DArrayOfFrequencyMagnitudePairs(x2, numberOfWindows)
	y2 = list(map(mappingFunction, xArray2))
	#print(len(y2))
	#print(y2)

	yCorrelated = np.correlate(y1, y2, mode='same')
	#print(yCorrelated)
	indices = np.argmax(yCorrelated)
	#print(indices)

	return indices


#x is signal in time domain, scipy.fft
# def bandpassFilter(x, f0, f1):
	
# 	#FFT
# 	X = scipy.fft(x)

# 	#Filter


# 	#IFFT


def butter_bandpass(lowcut, highcut, fs, order=5):
	nyq = 0.5 * fs
	low = lowcut / nyq
	high = highcut / nyq
	b, a = butter(order, [low, high], btype='band')
	return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
	b, a = butter_bandpass(lowcut, highcut, fs, order=order)
	y = lfilter(b, a, data)
	return y

def convertWindowIndexIntoSeconds(windowIndex):
	return windowIndex * (1.0/kSamplingRate) * (kWindowLength) * (kWindowShift/kWindowLength)

startEndArray = [(72324, 138915), 
					(101871, 156555), 
					(105840, 167580), 
					(19404, 102753),
					(315315, 403515), 
					(295470, 357210), 
					(275625, 379260), 
					(110250, 178605),
					(330750, 412335), 
					(168903, 244755),
				]

if __name__ == '__main__':

	arrayOfSetsOfPeaks = []
	#arrayOfSetsOfFilteredPeaks = []
	arrayOfCountOfWindows = []
	#arrayOf10SecClips = []

	for i in range(1,11):
		rate, data = scipy.io.wavfile.read('Data/clips/' + str(i) + '.wav')
	# Strip out the stereo channel if present
		if (len(data.shape) > 1):
			data = data[:,0]


	# Get just the first 10 seconds as our audio signal
		x = data[0:10*rate]
		#arrayOf10SecClips.append(x)

		#xFiltered = butter_bandpass_filter(x, kLowCut, kHighCut, kSamplingRate)


		X = stft(x)
		#XFiltered = stft(xFiltered)
		#print(X.shape)
		plot_transform(X)
	# Save the figure we just plotted as a .png
		pylab.savefig('spectrogram' + str(i) + '.png')
		#plot_transform(XFiltered)
	# Save the figure we just plotted as a .png
		#pylab.savefig('filteredSpectrogram' + str(i) + '.png')
		
	# Plot the peaks
		numberOfWindows = X.shape[0]
		arrayOfCountOfWindows.append(numberOfWindows)
		#print (numberOfWindows)
		s = setOfPeaksForTimeFrequencyDataWithMagnitude(X)
		arrayOfSetsOfPeaks.append(s)
		
		# sFiltered = setOfPeaksForTimeFrequencyDataWithMagnitude(XFiltered)
		# arrayOfSetsOfFilteredPeaks.append(sFiltered)

		#print (s)
		#peaks = list(map(lambda x: (x[0], x[1]), s))
		#print (peaks)
		plot_peaks(s, startEndArray[i-1][0]/kWindowLength * kWindowLength/kWindowShift,startEndArray[i-1][1]/kWindowLength * kWindowLength/kWindowShift)
		pylab.savefig('peaks' + str(i) + '.png')

		#print (peaks)
		# plot_peaks(sFiltered, startEndArray[i-1][0]/kWindowLength * kWindowLength/kWindowShift,startEndArray[i-1][1]/kWindowLength * kWindowLength/kWindowShift)
		# pylab.savefig('filteredPeaks' + str(i) + '.png')

		#yLargestMagnitude = map(largestMagnitude, s)

	


	clipList = range(1,5)
	print("**********************************")
	print("#LOWESTFREQUENCYINDEX #NOFILTER\n")
	#returns how many windows in the future 2 should be scaled
	for (i1, i2) in itertools.combinations(clipList, 2):
		windowOffset = crossCorrelate(arrayOfSetsOfPeaks[i1-1], arrayOfSetsOfPeaks[i2-1], arrayOfCountOfWindows[i1], lowestFrequencyIndex)-(arrayOfCountOfWindows[i1]/2)
		seconds = convertWindowIndexIntoSeconds(windowOffset)
		print("To Align {:d} with {:d}, shift {:d} by {:f} seconds.".format(i1, i2, i1, seconds))
	print("**********************************\n")

	print("**********************************")
	print("#HIGHESTFREQUENCYINDEX #NOFILTER\n")
	#returns how many windows in the future 2 should be scaled
	for (i1, i2) in itertools.combinations(clipList, 2):
		windowOffset = crossCorrelate(arrayOfSetsOfPeaks[i1-1], arrayOfSetsOfPeaks[i2-1], arrayOfCountOfWindows[i1], highestFrequencyIndex)-(arrayOfCountOfWindows[i1]/2)
		seconds = convertWindowIndexIntoSeconds(windowOffset)
		print("To Align {:d} with {:d}, shift {:d} by {:f} seconds.".format(i1, i2, i1, seconds))
	print("**********************************\n")

	print("**********************************")
	print("#LARGESTMAGNITUDE #NOFILTER\n")
	#returns how many windows in the future 2 should be scaled
	for (i1, i2) in itertools.combinations(clipList, 2):
		windowOffset = crossCorrelate(arrayOfSetsOfPeaks[i1-1], arrayOfSetsOfPeaks[i2-1], arrayOfCountOfWindows[i1], largestMagnitude)-(arrayOfCountOfWindows[i1]/2)
		seconds = convertWindowIndexIntoSeconds(windowOffset)
		print("To Align {:d} with {:d}, shift {:d} by {:f} seconds.".format(i1, i2, i1, seconds))
	print("**********************************\n")

	print("**********************************")
	print("#SUMOFMAGNITUDES #NOFILTER\n")
	#returns how many windows in the future 2 should be scaled
	for (i1, i2) in itertools.combinations(clipList, 2):
		windowOffset = crossCorrelate(arrayOfSetsOfPeaks[i1-1], arrayOfSetsOfPeaks[i2-1], arrayOfCountOfWindows[i1], sumOfMagnitudes)-(arrayOfCountOfWindows[i1]/2)
		seconds = convertWindowIndexIntoSeconds(windowOffset)
		print("To Align {:d} with {:d}, shift {:d} by {:f} seconds.".format(i1, i2, i1, seconds))
	print("**********************************\n")





	# print("**********************************")
	# print("#LOWESTFREQUENCYINDEX #FILTERED\n")
	# #returns how many windows in the future 2 should be scaled
	# for (i1, i2) in itertools.combinations(clipList, 2):
	# 	windowOffset = crossCorrelate(arrayOfSetsOfFilteredPeaks[i1-1], arrayOfSetsOfFilteredPeaks[i2-1], arrayOfCountOfWindows[i1], lowestFrequencyIndex)-(arrayOfCountOfWindows[i1]/2)
	# 	seconds = convertWindowIndexIntoSeconds(windowOffset)
	# 	print("To Align {:d} with {:d}, shift {:d} by {:f} seconds.".format(i1, i2, i1, seconds))
	# print("**********************************\n")

	# print("**********************************")
	# print("#HIGHESTFREQUENCYINDEX #FILTERED\n")
	# #returns how many windows in the future 2 should be scaled
	# for (i1, i2) in itertools.combinations(clipList, 2):
	# 	windowOffset = crossCorrelate(arrayOfSetsOfFilteredPeaks[i1-1], arrayOfSetsOfFilteredPeaks[i2-1], arrayOfCountOfWindows[i1], highestFrequencyIndex)-(arrayOfCountOfWindows[i1]/2)
	# 	seconds = convertWindowIndexIntoSeconds(windowOffset)
	# 	print("To Align {:d} with {:d}, shift {:d} by {:f} seconds.".format(i1, i2, i1, seconds))
	# print("**********************************\n")

	# print("**********************************")
	# print("#LARGESTMAGNITUDE #FILTERED\n")
	# #returns how many windows in the future 2 should be scaled
	# for (i1, i2) in itertools.combinations(clipList, 2):
	# 	windowOffset = crossCorrelate(arrayOfSetsOfFilteredPeaks[i1-1], arrayOfSetsOfFilteredPeaks[i2-1], arrayOfCountOfWindows[i1], largestMagnitude)-(arrayOfCountOfWindows[i1]/2)
	# 	seconds = convertWindowIndexIntoSeconds(windowOffset)
	# 	print("To Align {:d} with {:d}, shift {:d} by {:f} seconds.".format(i1, i2, i1, seconds))
	# print("**********************************\n")

	# print("**********************************")
	# print("#SUMOFMAGNITUDES #FILTERED\n")
	# #returns how many windows in the future 2 should be scaled
	# for (i1, i2) in itertools.combinations(clipList, 2):
	# 	windowOffset = crossCorrelate(arrayOfSetsOfFilteredPeaks[i1-1], arrayOfSetsOfFilteredPeaks[i2-1], arrayOfCountOfWindows[i1], sumOfMagnitudes)-(arrayOfCountOfWindows[i1]/2)
	# 	seconds = convertWindowIndexIntoSeconds(windowOffset)
	# 	print("To Align {:d} with {:d}, shift {:d} by {:f} seconds.".format(i1, i2, i1, seconds))
	# print("**********************************\n")


	# Wait for the user to continue (exiting the script closes all figures)
	input('Press [Enter] to finish')