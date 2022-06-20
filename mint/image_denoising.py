"""Functions used to reduce noise on individual video frames.

    tophat applies a white top-hat transform to each individual frame.
    wavelet_denoising applies a wavelet filter.
    lowpass1 and lowpass2 apply a low pass filter to the frames.
    wavelet_denoising applies a wavelet filter to the frames.
"""
import cv2
import numpy as np
from scipy import signal


def tophat(parameters,processed_frames,i,frames):
    """Applies top-hat transform to each individual frames. Inputs a dictionary, a NumPy array, an integer, and a NumPy array. Returns a NumPy array.

                Uses cv2 MORPH_TOPHAT function to remove artifacts.
                parameters is a dictionary containing calculation parameters. This function uses the 'separation' key. Refer to script.py for its use.
                processed_frames is an empty NumPy array, shaped according to the number of frames, rows and columns of the current file.
                i is the index of the for loop processing each frame.
                frames is a NumPy array containing frames from the current file.
                """              
    kernelC = cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(parameters['separation'],parameters['separation']))
    processed_frames[i] = cv2.morphologyEx(frames[i],cv2.MORPH_TOPHAT,kernelC)

    return processed_frames

def lowpass1():
    """Defines an array for low pass filtering, returns a NumPy array.
    """
    h0 = 3./8
    h1 = 1./4
    h2 = 1./16
    lowpass = np.array([h2,h1,h0,h1,h2])
    lowpassV = lowpass.reshape(-1,1)
    lowpassH = lowpass.reshape(1,-1)
    array1 = np.dot(lowpassV,lowpassH)
    return array1

def lowpass2():
    """Defines an array for low pass filtering, returns a NumPy array.
    """
    h0 = 3./8
    h1 = 1./4
    h2 = 1./16
    lowpass = np.array([h2,0,h1,0,h0,0,h1,0,h2])
    lowpassV = lowpass.reshape(-1,1)
    lowpassH = lowpass.reshape(1,-1)
    array2 = np.dot(lowpassV,lowpassH)
    return array2

def wavelet(array,array1,array2):
    """Applies a wavelet filter. Inputs three NumPy arrays, returns a NumPy array.
    """
    lc1 = signal.convolve2d(array,array1,mode='same')
    lc2 = signal.convolve2d(lc1,array2,mode='same')
    new_array = lc1 - lc2
    return new_array

def wavelet_denoising(processed_frames,i):
    """Applies wavelet denoising to each individual frames. Inputs a NumPy array and an integer, returns a NumPy array.

        Uses two low-pass filers and scipy.signal.convolve2d to remove background noise.
        processed_frames is a NumPy array, either filled with raw frames or frames that went through top-hat transform.
        i is the index of the for loop processing each frame.
        """ 
    lp1 = lowpass1()
    lp2 = lowpass2()

    processed_frames[i] = wavelet(processed_frames[i],lp1,lp2)

    return processed_frames

def line_average(frames):
    """Performs line-by-line frame average. Inputs and returns an array of frames.
        Each line of even frames is averaged with the corresponding line of the following odd frames.
    """
    avg_frames = np.empty((int(frames.shape[0]/2),frames.shape[1],frames.shape[2]))

    for i in range(0,(frames.shape[0]-2),2):
        for j in range(frames.shape[1]):
            avg_frames[(int(i/2)),j,] = ((frames[i,j,]+frames[i+1,j,])/2)

    return avg_frames

def frame_average(frames):
    """Performs line-by-line frame average. Inputs and returns an array of frames.
        Each line of even frames is averaged with the corresponding line of the following odd frames.
    """
    avg_frames = np.empty((int(frames.shape[0]/2),frames.shape[1],frames.shape[2]))

    for i in range(0,(frames.shape[0]-2),2):
        avg_frames[(int(i/2))] = ((frames[i]+frames[i+1])/2)

    return avg_frames

def frame_accu(frames):
    """Performs line-by-line frame average. Inputs and returns an array of frames.
        Each line of even frames is averaged with the corresponding line of the following odd frames.
    """
    avg_frames = np.empty((int(frames.shape[0]/2),frames.shape[1],frames.shape[2]))

    for i in range(0,(frames.shape[0]-2),2):
        avg_frames[(int(i/2))] = ((frames[i]+frames[i+1]))

    return avg_frames