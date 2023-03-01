"""Functions used to reduce noise on individual video frames."""

import cv2
import numpy as np
from scipy import signal

def tophat(separation,frame):
    """Applies top-hat transform.

    Frame-by-frame top-hat filtering with `cv2.MORPH_TOPHAT`. Removes artifacts.

    :param separation: Minimum distance (in pixels) between features.
    :type separation: dict
    :param frame: 2D array of unfiltered data.
    :type frame: NumPy array
    :return: 2D array of filtered data.
    :rtype: NumPy array
    """    

    kernelC = cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(separation,separation))
    processed_frame = cv2.morphologyEx(frame,cv2.MORPH_TOPHAT,kernelC)

    return processed_frame

def lowpass():
    """Defines a pair of arrays for low pass filtering.

    :return: Arrays for low pass filtering.
    :rtype: NumPy array
    """    

    h0 = 3./8
    h1 = 1./4
    h2 = 1./16

    lowpass1 = np.array([h2,h1,h0,h1,h2])
    lowpassV1 = lowpass1.reshape(-1,1)
    lowpassH1 = lowpass1.reshape(1,-1)
    array1 = np.dot(lowpassV1,lowpassH1)

    lowpass2 = np.array([h2,0,h1,0,h0,0,h1,0,h2])
    lowpassV2 = lowpass2.reshape(-1,1)
    lowpassH2 = lowpass2.reshape(1,-1)
    array2 = np.dot(lowpassV2,lowpassH2)

    return array1,array2

def wavelet(frame):
    """Applies wavelet denoising.

    Uses two low-pass filers and `scipy.signal.convolve2d` to remove background noise.

    :param frame: 2D array of unfiltered data.
    :type frame: NumPy array
    :return: 2D array of filtered data.
    :rtype: NumPy array
    """    

    lp1 = lowpass()[0]
    lp2 = lowpass()[1]

    lc1 = signal.convolve2d(frame,lp1,mode='same')
    lc2 = signal.convolve2d(lc1,lp2,mode='same')
    processed_frame = lc1 - lc2

    return processed_frame
