"""Functions used to reduce noise on individual video frames."""

import os
import cv2
import numpy as np
from scipy import signal
from joblib import Parallel, delayed


def tophat(separation: int, frame: np.ndarray) -> np.ndarray:
    """Applies top-hat transform.

    Top-hat filtering with `cv2.MORPH_TOPHAT`. Removes artifacts.

    :param separation: Minimum distance (in pixels) between features.
    :type separation: dict
    :param frame: 2D array of unfiltered data.
    :type frame: np.ndarray
    :return: 2D array of filtered data.
    :rtype: np.ndarray
    """

    kernelC = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (separation, separation))
    processed_frame = cv2.morphologyEx(frame, cv2.MORPH_TOPHAT, kernelC)

    return processed_frame

def lowpass() -> tuple[np.ndarray, np.ndarray]:
    """Defines a pair of arrays for low pass filtering.

    :return: Arrays for low pass filtering.
    :rtype: np.ndarray
    """

    h0 = 3./8
    h1 = 1./4
    h2 = 1./16

    lowpass1 = np.array([h2, h1, h0, h1, h2])
    lowpassV1 = lowpass1.reshape(-1, 1)
    lowpassH1 = lowpass1.reshape(1, -1)
    array1 = np.dot(lowpassV1, lowpassH1)

    lowpass2 = np.array([h2, 0, h1, 0, h0, 0, h1, 0, h2])
    lowpassV2 = lowpass2.reshape(-1, 1)
    lowpassH2 = lowpass2.reshape(1, -1)
    array2 = np.dot(lowpassV2, lowpassH2)

    return array1, array2

def wavelet(frame: np.ndarray) -> np.ndarray:
    """Applies wavelet denoising.

    Uses two lowpass filers and `scipy.signal.convolve2d` to remove background noise.

    :param frame: 2D array of unfiltered data.
    :type frame: np.ndarray
    :return: 2D array of filtered data.
    :rtype: np.ndarray
    """

    lp1 = lowpass()[0]
    lp2 = lowpass()[1]

    lc1 = signal.convolve2d(frame, lp1, mode='same')
    lc2 = signal.convolve2d(lc1, lp2, mode='same')
    processed_frame = lc1 - lc2

    return processed_frame

def filtering(frames: np.ndarray, settings: dict, parameters: dict) -> np.ndarray:
    """Performs frame-by-frame filtering.

    :param frames: 3D array of video frames.
    :type frames: np.ndarray
    :param settings: Dictionary of filtering settings.
    :type settings: dict
    :param parameters: Dictionary of filtering parameters.
    :type parameters: dict
    :return: 3D array of filtered video frames.
    :rtype: np.ndarray
    """

    for i in range(len(frames)):
        if settings['tophat']:
            frames[i] = tophat(parameters['separation'], frames[i])

        if settings['wavelet']:
            frames[i] = wavelet(frames[i])

    return frames

def filtering_p(frames: np.ndarray, settings: dict, parameters: dict) -> np.ndarray:
    """Performs frame-by-frame filtering using parallel processing.

    :param frames: 3D array of video frames.
    :type frames: np.ndarray
    :param settings: Dictionary of filtering settings.
    :type settings: dict
    :param parameters: Dictionary of filtering parameters.
    :type parameters: dict
    :return: 3D array of filtered video frames.
    :rtype: np.ndarray
    """

    if settings['tophat']:

        tophat_gen = Parallel(n_jobs=os.cpu_count(),
                              return_as='generator')(delayed(tophat)
                                                    (parameters['separation'], frame)
                                                    for frame in frames)

        for i, frame in zip(range(len(frames)), tophat_gen):
            frames[i] = frame

    if settings['wavelet']:

        wavelet_gen = Parallel(n_jobs=os.cpu_count(),
                               return_as='generator')(delayed(wavelet)(frame)
                                                       for frame in frames)

        for i, frame in zip(range(len(frames)), wavelet_gen):
            frames[i] = frame

    return frames