"""Functions used to reduce noise on individual video frames."""

import os
import cv2
import numpy as np
from scipy import signal
from joblib import Parallel, delayed

def tophat(
    separation: int,
    frame: np.ndarray,
    ) -> np.ndarray:
    """
    Applies top-hat transform to input image.

    Top-hat filtering with `cv2.MORPH_TOPHAT`. Removes artifacts.

    Parameters
    ----------
    separation : int
        Minimum distance (in pixels) between features.
    frame : np.ndarray
        2D array of unfiltered data.

    Returns
    -------
    np.ndarray
        2D array of filtered data.
    """

    kernelC = np.ones((separation, separation), dtype=int)
    return frame - cv2.erode(frame, kernelC)

def lowpass(
    ) -> tuple[np.ndarray, np.ndarray]:
    """
    Defines a pair of arrays for low pass filtering.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Arrays for low pass filtering.
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

def wavelet(
        frame: np.ndarray,
        ) -> np.ndarray:
    """
    Applies wavelet denoising to a single frame.

    Uses two lowpass filters and `scipy.signal.convolve2d` to remove background noise.

    Parameters
    ----------
    frame : np.ndarray
        2D array of unfiltered data.

    Returns
    -------
    np.ndarray
        2D array of filtered data.
    """

    lp1 = lowpass()[0]
    lp2 = lowpass()[1]

    lc1 = signal.convolve2d(frame, lp1, mode='same')
    lc2 = signal.convolve2d(lc1, lp2, mode='same')

    return lc1 - lc2

def filtering(
        frames: np.ndarray,
        settings: dict,
        parameters: dict,
        ) -> np.ndarray:
    """
    Perform frame-by-frame filtering.

    Parameters
    ----------
    frames : np.ndarray
        3D array of video frames.
    settings : dict
        Dictionary of filtering settings.
    parameters : dict
        Dictionary of filtering parameters.

    Returns
    -------
    np.ndarray
        3D array of filtered video frames.
    """

    if len(frames.shape) == 2:
        frames = np.expand_dims(frames, axis=0)

    for i, frame in enumerate(frames):
        if settings['tophat']:
            frames[i] = tophat(parameters['separation'], frame)

        if settings['wavelet']:
            frames[i] = wavelet(frame)

    return frames

def array_filtering(
        frames: np.ndarray,
        settings: dict,
        parameters: dict,
        ) -> np.ndarray:
    """
    Perform frame-by-frame filtering.

    Parameters
    ----------
    frames : np.ndarray
        3D array of video frames.
    settings : dict
        Dictionary of filtering settings.
    parameters : dict
        Dictionary of filtering parameters.

    Returns
    -------
    np.ndarray
        3D array of filtered video frames.
    """

    if settings['tophat']:
        frames = np.array([tophat(separation=parameters['separation'], frame=frame)
                           for frame in frames])

    if settings['wavelet']:
        frames = np.array([wavelet(frame=frame) for frame in frames])

    return frames

def filtering_p(
        frames: np.ndarray,
        settings: dict,
        parameters: dict,
        ) -> np.ndarray:
    """
    Perform frame-by-frame filtering.

    Parameters
    ----------
    frames : np.ndarray
        3D array of video frames.
    settings : dict
        Dictionary of filtering settings.
    parameters : dict
        Dictionary of filtering parameters.

    Returns
    -------
    np.ndarray
        3D array of filtered video frames.
    """

    if settings['tophat']:

        tophat_gen = Parallel(n_jobs=os.cpu_count(),
                              return_as='generator')(delayed(tophat)
                                                    (parameters['separation'], frame)
                                                    for frame in frames)
        for i, frame in enumerate(tophat_gen):
            frames[i] = frame

    if settings['wavelet']:

        wavelet_gen = Parallel(n_jobs=os.cpu_count(),
                               return_as='generator')(delayed(wavelet)(frame)
                                                       for frame in frames)
        for i, frame in enumerate(wavelet_gen):
            frames[i] = frame

    return frames