import scipy
import scipy.io.wavfile
import pylab
import matplotlib
import numpy as np

kWindowSize = 20
kWindowLength = 4096
kWindowShift = 2048

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

def setOfPeaksForTimeFrequencyData(X):
	setOfMaxima = set()
	Xdims = X.shape
	for i in range(0, Xdims[0]):
		for j in range(0, Xdims[1]):
			if isMaxima(i, j, X):
				setOfMaxima.add((i,j))
	return setOfMaxima


# Plot a list of peaks in the form [(s1, f1), (s2, f2), ...]
def plot_peaks(peak_list,start,end):
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(1,1,1)    

    filteredPeakList = list(filter(lambda x: (x[0]>=start and x[0]<=end), peak_list))

    #s_list, f_list = zip(*filteredPeakList)

    s_list, f_list = zip(*peak_list)
    #print (s_list)



    matplotlib.pyplot.plot(s_list, f_list, 'bo')    
    ymin, ymax = ax.get_ylim()    
    ax.vlines((start,end),ymin,ymax,'red')
    matplotlib.pyplot.xlabel('Window index')
    matplotlib.pyplot.ylabel('Transform coefficient')
  

startEndArray = [(72324, 138915), (101871, 156555), (105840, 167580), (19404, 102753)]

if __name__ == '__main__':
	for i in range(1,5):
		rate, data = scipy.io.wavfile.read('Data/clips/' + str(i) + '.wav')
	# Strip out the stereo channel if present
		if (len(data.shape) > 1):
			data = data[:,0]

	# Get just the first 10 seconds as our audio signal
		x = data[0:10*rate]

		X = stft(x)
		#print(X.shape)
		plot_transform(X)
	# Save the figure we just plotted as a .png
		pylab.savefig('spectrogram' + str(i) + '.png')
		
	# Plot the peaks
		s = setOfPeaksForTimeFrequencyData(X)
		plot_peaks(s, startEndArray[i-1][0]/kWindowLength * kWindowLength/kWindowShift,startEndArray[i-1][1]/kWindowLength * kWindowLength/kWindowShift)
		pylab.savefig('peaks' + str(i) + '.png')

	# Wait for the user to continue (exiting the script closes all figures)
	input('Press [Enter] to finish')