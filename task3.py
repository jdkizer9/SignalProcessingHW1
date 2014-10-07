import os
import scipy
import scipy.io.wavfile
import pylab
import matplotlib
import numpy as np


#Task 1 :
def haar(X):
    #convert items into float first otherwise calculations will be rounded
    X = list(map(float, X))
    #Turn each adjacent pair into an (avg,diff) pair
    avg_result = []
    diff_result = []
    for i in range(0,len(X)):
        # read every other x to extract pairs
        if i % 2 != 0: continue
        avg = ( X[i] + X[i+1] ) / 2
        diff = X[i] - avg
        avg_result.append(avg)
        diff_result.append(diff)
    
    #Recurse on the left half of the output
    if len(avg_result) > 1:
        avg_result = haar(avg_result)
    
    # averages on the left and differences on the right
    return avg_result + diff_result
    

# #haar([9,7,3,5]) should return [6,2,1,-1]
# print (haar([9,7,3,5]))

# #haar([1,2,3,4,5,6,7,8]) should return [4.5, -2.0, -1.0, -1.0, -0.5, -0.5, -0.5, -0.5]
# print (haar([1,2,3,4,5,6,7,8]))

# #haar([9,7,3,5,0,8,2,9]) should return [5.375, 0.625, 2.0, -0.75, 1.0, -1.0, -4.0, -3.5]
# print (haar([9,7,3,5,0,8,2,9]))


print ('inverse haar')
#inverse haar
def inverse_haar(X):
    #convert items into float first otherwise calculations will be rounded
    X = list(map(float, X))
    
    #source: http://people.sc.fsu.edu/~jburkardt/cpp_src/haar/haar.cpp
    #initialize the output array with same len as input
    y = list(range(len(X)))
    
    #initialize a counter that increments in *2
    k = 1
    while k * 2 <= len(X):
        for i in range(k):
            #form the first and second terms from the average and difference
            y[2*i] = X[i] +X[i+k]
            y[2*i+1] = X[i] -X[i+k]
        for i in range(k*2):
            #replace the terms in input X, next loop would require them in calculation
            X[i] = y[i]
        k= k*2
    
    return y

    
# #inverse_haar([6,2,1,-1]) should return [9,7,3,5]
# print (inverse_haar([6,2,1,-1]))

# #inverse_haar([4.5, -2.0, -1.0, -1.0, -0.5, -0.5, -0.5, -0.5]) should give me [1,2,3,4,5,6,7,8]
# print (inverse_haar([4.5, -2.0, -1.0, -1.0, -0.5, -0.5, -0.5, -0.5]))

# #inverse_haar([5.375, 0.625, 2.0, -0.75, 1.0, -1.0, -4.0, -3.5])  should give me [9,7,3,5,0,8,2,9]
# print (inverse_haar([5.375, 0.625, 2.0, -0.75, 1.0, -1.0, -4.0, -3.5]))    

#Task 2

def graph_accel(signal):
	pylab.ion()
	pylab.figure()
	pylab.plot(signal)
	pylab.ioff()

#trim the input to the largest 2^n-sized subsequence.
def trimSignal(x):
    largest_n =  np.floor(np.log2(len(x)))
    largest_len = int(2**largest_n)
    
    #this trims all the right-hand side signals
    return x[0:largest_len]
    
    #this trims all the left-hand side signals
    #return x[len(x)-largest_len::]
    

filedir = os.path.dirname(__file__) + '/Data/'


f = open(filedir + 'accel.csv', 'r')

# Get just the displacement in the x coordinate
original_x = [float(line.split(',')[6]) for line in f]


#A. Haar transformed data
#First, some basic applications of the Haar transform to this data.
#1. Use the function you wrote in Task 1 to transform the data. Recall that Haar requires input
#which has a length of 2^n, so you should first trim the input to the largest 2^n-sized
#subsequence.
#2. Graph the result using graph_accel().

#2 A1). Trim the signal to 2^n length and transform to haar wavelet
x = trimSignal(original_x)
graph_accel(x)
pylab.savefig(filedir + 'accel_trimmedsignal.png')

x = haar(x)

#2 A2). Graph the result using graph_accel().
# Graph it, and save figure as a .png
graph_accel(x)
pylab.savefig(filedir + 'accel_haar.png')

#B. Haar transforms and “edges”
#One particular reason for using the Haar transform over other wavelet transforms is that the basic
#operation it uses is taking pairwise differences. Pairwise differences will be large when the signal is
#changing rapidly, and small when the signal is constant, so we can use the Haar transform to find rapid
#changes in the signal, or “steps”.
#First, look at just the second half (i.e., the last 2n-1 entries) of the Haar-transformed signal. These
#correspond to pairwise differences of individual entries of the original signal. These are typically referred
#to as the first-order differences, so we'll denote this by X1.
#Write a function to find the 5 biggest positive entries of the pairwise differences X1, and the 5 most
#negative entries. Briefly (one or two sentences) describe what is happening to the original signal at the
#locations corresponding to these points.

#extract first order difference x_one from x
x_one = x[int(len(x)-len(x)/2)::]

#print x_one
#print len(x_one)

sorted_x_one = np.sort(x_one)

print ('bottom/smallest 5 x_one:', sorted_x_one[0:5])
print ('top/largest 5 x_one:', sorted_x_one[len(sorted_x_one)-5::])

#The largest and smallest values 0.2195 and -1.3355 corresponds to the abrupt drop and rise that looks like a square wave at around row 112 to 149
#The other peaks are caused by abrupt changes around row 1, 30 and 211

#C. Smoothing
#Loosely speaking, the second half of the Haar transformed signal captures rapid changes in the original
#signal, while the earlier entries capture more slowly varying changes. This suggests that we can get a
#smoothed version of the input by ignoring the information in the second half of the Haar-transformed
#signal.
#In particular, we already have a way of inverting the Haar-transform, the function inverse_haar(X) that
#you wrote earlier. So, we can just throw away the information corresponding to rapid changes, and take
#the inverse transform of the result. Let's see what happens when we do this.
#1. First, take the original signal x and compute its Haar-transform X. Call inverse_haar(X) and
#graph the result. Compare the result to the original signal (hint: they should be identical).
#2. Start with the Haar-transformed signal X, and truncate X to the first 2^(n-1) entries. Take the
#inverse Haar-transform of this truncated X, and graph the result.
#3. Do the same thing with the first 2^(n-k) entries, for k = 2, 3 and 4, graphing each result. Briefly
#describe (1 to 2 sentences) what happens as k increases.

#1.verify that inverse_haar function produces the same result
trimmed_x = trimSignal(original_x)
transformed_x = haar(trimmed_x)
inverse_transformed_x = inverse_haar(transformed_x)
# The sum of differences between original signal and inverse signal should be very close to 0, and therefore numerically equivalent 
# These differences arise from using float types for calculation, which is faster but sacrifices accuracy. For non-scientific calculation use, this error is acceptable.
print ('Sum of differences between original signal and inverse signal: ', np.sum(np.absolute(np.subtract(trimmed_x,inverse_transformed_x))))

#2. Start with the Haar-transformed signal X, and truncate X to the first 2^(n-1) entries. Take the
#inverse Haar-transform of this truncated X, and graph the result.
truncated_x = transformed_x[0:int(len(transformed_x)/2)]
graph_accel(inverse_haar(truncated_x))
pylab.savefig(filedir + 'accel_haar_firsthalf_k_1.png')

#3. Do the same thing with the first 2^(n-k) entries, for k = 2, 3 and 4, graphing each result. Briefly
#describe (1 to 2 sentences) what happens as k increases.
for k in [2,3,4]:
    truncated_x = transformed_x[0:int(len(transformed_x)/(2**k))]
    graph_accel(inverse_haar(truncated_x))
    pylab.savefig(filedir + 'accel_haar_firsthalf_k_{}.png'.format(k))

#Looking at the different charts it appears as k increases, the chart loses details. With fewer and fewer data points trying to represent the entire chart, only the most exaggerated features such as large drops and rises are apparent.


# Task 3:

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

		s = setOfPeaksForTimeFrequencyData(X)

	# Plot some dummy peaks
		#plot_peaks([(100, 50), (200, 87), (345, 20)],150,200)
		#s = set()
		#s.add((100, 50))
		#s.add((200, 87))
		#s.add((345, 20))
		#print(s)
		plot_peaks(s, startEndArray[i-1][0]/kWindowLength * kWindowLength/kWindowShift,startEndArray[i-1][1]/kWindowLength * kWindowLength/kWindowShift)
		pylab.savefig('peaks' + str(i) + '.png')

	# Wait for the user to continue (exiting the script closes all figures)
	input('Press [Enter] to finish')