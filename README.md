Tracking
===

Tracking is a Python script used to extract intraneuronal transport parameters from microscopy data.

It relies on [**TrackPy**](https://github.com/soft-matter/trackpy) to generate trajectories from video files, from which transport parameters are then extracted and statistically tested.

Its purpose is to automate workflow to the point where an input of raw video files results in an output of p-values and graphs.

Please refer to the [**documentation**]() for further information.



Quickstart
---

### **CLI**
Download the lastest release and extract its contents to a dedicated folder.

Edit calculation parameters and settings as needed in `script.py`.

Run `script.py` through command line or any IDE.


### **GUI** (work in progress)
Download the latest release and extract its contents to a dedicated folder.

Run `gui.py`.