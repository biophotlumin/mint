"""Functions used to reduce noise on individual video frames.
"""

import cv2
import numpy as np
from scipy import signal

def tophat(parameters,processed_frame):
    """Applies top-hat transform.

    Frame-by-frame top-hat filtering with `cv2.MORPH_TOPHAT`. Removes artifacts.

    :param parameters: Dictionary containing the minimum distance (in pixels) between features under the `'separation'` key.
    :type parameters: dict
    :param processed_frames: Empty 3D array of identical shape to the raw frames array.
    :type processed_frames: NumPy array
    :param i: Index of the current frame.
    :type i: int
    :param frames: 3D array of raw frames over time.
    :type frames: NumPy array
    :return: 3D array of processed frames over time.
    :rtype: NumPy array
    """    

    kernelC = cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(parameters['separation'],parameters['separation']))
    processed_frame = cv2.morphologyEx(processed_frame,cv2.MORPH_TOPHAT,kernelC)

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

def wavelet(processed_frame):
    """Applies wavelet denoising.

    Uses two low-pass filers and `scipy.signal.convolve2d` to remove background noise.

    :param processed_frames: 3D array of frames over time.
    :type processed_frames: NumPy array
    :param i: Index of the current frame.
    :type i: int
    :return: 3D array of processed frames over time.
    :rtype: NumPy array
    """    

    lp1 = lowpass()[0]
    lp2 = lowpass()[1]

    lc1 = signal.convolve2d(processed_frame,lp1,mode='same')
    lc2 = signal.convolve2d(lc1,lp2,mode='same')
    processed_frame = lc1 - lc2

    return processed_frame
