"""
Example of Fourier Spectrum Analysis using FFT. 

Generates a sin wave in the time domain and uses FFT to 
convert signal to frequency domain using numpy fft. 
"""


import numpy as np
from math import pi
import matplotlib.pyplot as plt

# generate sin wave as sample data in time domain (amplitude vs time)
pts = 128
sinFreq = 5.0
t = np.linspace(0,1,pts)
amp = np.sin(2*pi*sinFreq*t)

# do FFT on time domain signal to get freq domain (freq spectrum vs freq)
s = abs(np.fft.fftn(amp))   # full two-sided spectrum, take abs since FFT gives complex values
spectrum = s[range(pts/2)]  # take only one side of spectrum
fmin = min(spectrum)
fmax = max(spectrum)
npts = len(spectrum)
freq = np.linspace(fmin,fmax,npts)

# plot time domain signal
plt.subplot(211)
plt.title('FFT of Time Signal to Freq Signal')
plt.xlabel('Time (sec)')
plt.ylabel('Amplitude')
plt.plot(t,amp)

# plot freq domain signal
plt.subplot(212)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Spectrum')
plt.plot(freq,spectrum)