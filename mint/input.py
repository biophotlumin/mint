#Imports

import abc
import warnings
from pathlib import Path

import nd2
import imageio
import numpy as np

from .utils import Path_type

try:
    import imagej
    from jpype.types import JClass

    IJ_INSTALLED = True

except ImportError:

    IJ_INSTALLED = False

class BaseReader(abc.ABC):
    """
    Abstract base class for video readers.

    All derived classes should implement the `get_frames` method to return
    the frames of the video as a 3D numpy array.
    """

    def __init__(self):
        self.name = 'Base reader'

    @abc.abstractmethod
    def get_frames(self,
                   file_path: Path_type,
                   ) -> np.ndarray:
        """
        Extract the frames from the video.

        Parameters
        ----------
        file_path : str or pathlib.Path
            File path to the video.

        Returns
        -------
        frames : np.ndarray
            3D array of shape (number of frames, height, width) containing
            the frames of the video.
        """
        pass

    def get_shape(self,
                  frames: np.ndarray,
                  ) -> tuple[int, ...]:
        """
        Returns the shape of the frames.

        Parameters
        ----------
        frames : np.ndarray
            3D array of frames.

        Returns
        -------
        shape : tuple
            Shape of the frames.
        """
        return frames.shape

    def check_shape(self,
                    shape: tuple[int, ...],
                    dims: int = 3,
                    ) -> bool:
        """
        Check that the shape of the frames has the expected number of
        dimensions.

        Parameters
        ----------
        shape : tuple
            Shape of the frames.
        dims : int, optional
            Expected number of dimensions, by default 3.

        Returns
        -------
        bool
            True if the frames have the expected number of dimensions, False
            otherwise.
        """
        if len(shape) == dims:
            return True
        else:
            warnings.warn(f'Number of dimension(s) ({len(shape)})'
                          f'does not match expected {dims}')
            return False

    def return_frames(self,
                      file_path: Path_type,
                      ) -> np.ndarray:
        """
        Check and extract frames from the video.

        Parameters
        ----------
        file_path : str or pathlib.Path
            File path to the video.

        Returns
        -------
        frames : np.ndarray
            Returns the frames of the video.
        """

        try:
            frames = self.get_frames(file_path)
        except FileNotFoundError:
            warnings.warn(f'File {file_path} not found')
            return np.empty(0)

        if isinstance(frames, np.ndarray):
            if self.check_shape(self.get_shape(frames)) is True:
                return frames
            else:
                warnings.warn(f'Reading {Path(file_path).name} with {self.name}'
                              f' did not return an array of correct shape')
                return np.empty(0)
        else:
            warnings.warn(f'Reading {Path(file_path).name} with {self.name}'
                          f' did not return an array')
            return np.empty(0)

class ImageIOReader(BaseReader):

    def __init__(self):
        self.name = 'ImageIO reader'

    def get_frames(self,
                   file_path: Path_type,
                   ) -> np.ndarray:
        return imageio.volread(file_path)

class ND2Reader(BaseReader):

    def __init__(self):
        self.name = 'ND2 reader'

    def get_frames(self,
                   file_path: Path_type,
                   ) -> np.ndarray:
        return nd2.imread(file_path)

class BioFormatsReader(BaseReader):

    def __init__(self):
        self.name = 'BioFormats reader'
        self.jvm_started = False

    def start_JVM(self) -> JClass:
        self.jvm_started = True
        return imagej.init()

    def get_frames(self,
                   file_path: Path_type,
                   ) -> np.ndarray:

        file_path = str(file_path)

        if self.jvm_started is False:
            ij = self.start_JVM()
        jimage = ij.io().open(file_path)
        frames = ij.py.from_java(jimage)

        try:
            frames.shape
        except AttributeError:
            warnings.warn(f'Extension {file_path.split(".")[-1]} is not supported')

        return frames

def get_frames(
        file_path: Path_type,
        ) -> np.ndarray:
    """
    Loads frames from file.

    Parameters
    ----------
    file_path : Path_type
        File path.

    Returns
    -------
    frames : ndarray
        3D array of frames.
    """

    extension = Path(file_path).suffix

    if extension == '.tif':
        reader = ImageIOReader()
    elif extension == '.nd2':
        reader = ND2Reader()
    elif IJ_INSTALLED is True:
        reader = BioFormatsReader()
    else:
        warnings.warn(f'Extension {extension} is not supported')

    return reader.return_frames(file_path)

## TODO Assert float64 ?