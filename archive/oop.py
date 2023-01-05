import pandas as pd
import numpy as np

class Trajectory():

    def __init__(self) -> None:
        
        self.x = np.array()
        self.y = np.array()
        self.mass = np.array()
        self.frame = np.array()
        self.particle = np.array()
