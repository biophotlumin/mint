#Imports

import abc
import warnings
from pathlib import Path

import imageio
import nd2
import numpy as np

from .utils import Path_type

try:
    import imagej

    IJ_INSTALLED = True

except ImportError:

    IJ_INSTALLED = False

class BaseReader(abc.ABC):
    """
    Abstract base class for video readers.
    """

    def __init__(self):
        self.name = 'Base reader'

    @abc.abstractmethod
    def get_frames(self, file_path: Path_type) -> None:
        pass

    def get_shape(self, frames: np.ndarray) -> tuple:
        return frames.shape

    def check_shape(self, shape: tuple, dims: int = 3) -> bool:
        if len(shape) == dims:
            return True
        else:
            warnings.warn(f'Number of dimension(s) ({len(shape)})'
                          f'does not match expected {dims}')
            return False

    def return_frames(self, file_path: Path_type) -> np.ndarray:
        frames = self.get_frames(file_path)
        if isinstance(frames, np.ndarray):
            if self.check_shape(self.get_shape(frames)) is True:
                return frames
        else:
            warnings.warn(f'Reader {self.name} did not return an array')

class ImageIOReader(BaseReader):

    def __init__(self):
        self.name = 'ImageIO reader'

    def get_frames(self, file_path: Path_type) -> np.ndarray:
        return imageio.volread(file_path)

class ND2Reader(BaseReader):

    def __init__(self):
        self.name = 'ND2 reader'

    def get_frames(self, file_path: Path_type) -> np.ndarray:
        return nd2.imread(file_path)

class BioFormatsReader(BaseReader):

    def __init__(self):
        self.name = 'BioFormats reader'
        self.jvm_started = False

    def start_JVM(self):
        self.jvm_started = True
        return imagej.init()

    def get_frames(self, file_path: Path_type) -> np.ndarray:

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

def get_frames(file_path: Path_type) -> np.ndarray:
    """Get frames from a file.
    Reader selected based on file extension.

    :param file_path: File path
    :type file_path: Path_type
    :return: 3D array of frames
    :rtype: np.ndarray
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