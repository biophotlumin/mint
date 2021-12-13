********************
**User guide**
********************
.. contents:: :backlinks: None
 

**General information**
===============================

Input data must be in the form of .tif files, and should be placed in a dedicated folder. Proper folder structure is detailed in :ref:`Statistical analysis`.

Default parameters are provided for both CLI and GUI, but as they are very specific to our experimental setup please expect a lot of trial and errors 
before obtaining satisfactory results. 

You can use the "Test locate" function to try feature findind parameters.


**Parameters and settings**
===============================

**Preprocessing**
---------------------

Two preprocessing filters can be applied to the videos before going through trackpy :


* **Wavelet** denoising reduces background noise and improves signal-to-noise ratio.


* **Top-hat** transform (white) removes large spots and artifacts.

|

**Feature finding and linking**
-----------------------------------

This part is handled by `\trackpy <https://github.com/soft-matter/trackpy>`_. 
Please refer to their own documentation for more details on each parameter.


* 
  **Feature finding** :


  * **diameter** : Size of the features to be found. Must be odd.
  * **minmass** : Minimum integrated brightness.
  * **separation** : Minimum separation between features.

* 
  **Linking** : 


  * **search_range** : Maximum distance features can move between frames.
  * **memory** (optional) : Maximum number of frames a feature can disappear for and still be considered the same spot.
  * **adaptive_stop** (optional) : Bottom line search_range that trackpy will gradually moves towards if the original parameter results in an unsolvable calculation.
  * **adaptive_step** (optional) : Factor by which search_range will decrease until it reached adaptive_stop.
  * **stub_filtering** (setting, optional) : Filters trajectories based on length.

    * **stub_filtering** (parameter) : Minimum of points a trajectory must have to pass filtering.

* 
  **Additionally** : 


  * **MSD** (optional) : Computes the Mean Square Displacement (MSD) for each trajectory, and filters them accordingly.

    * **threshold** : MSD value below which a trajectory is discarded.

  * **rejoining** (optional) : Rejoins trajectories that weren't linked by trackpy in the first place. This function is kept deliberately stringent (i.e. each trajectory cannot be rejoined more than once) to avoid aberrant trajectories and false positives.

    * **threshold_t** : Maximum number of frames between the first and last points of two trajectories for them to be considered for rejoining.
    * **threshold_r** : Maximum distance (in pixels) between the first and last points of two trajectories for them to be considered for rejoining.

  * **SNR_estimation** : Calcualtes the signal-to-noise ratio.

    * **base_level** : Base noise level.

|

**Data extraction**
-----------------------

This part of the script calculates transport parameters from extracted trajectories.


* **r_conf_cut** : Cutoff level for confinement ratio calculation.
* **px** : Size of a pixel in µm.
* **dt** : Sampling rate, in seconds, as in the amount of time between two frames.
* **min_theoretical_precision** : Minimum theoretical precision, in nm.
* **sliding_window** : Size, in frames, of the sliding window along which the confinement ration is calculated.
* **polynomial_fit** (setting, optional) : Filters trajectories based on how well they fit to a 3rd degree polynom.

  * **threshold_poly3** : Tolerated deviation from the 3rd degree polynom.
  * **len_cutoff** : Size, in points, below which trajectories are eliminated.

* **minimization** (setting, optional) : Experimental function of convex minimization of acceleration used to smooth trajectories with a lot of spatial noise.

  * **sigma** : Estimated noise, in nm.

* **antero_retro** : Separates transport parameters into anterograde and retrograde categories. Note : this highly dependent on our original experimental setup and will most likely not work elsewhere. 

|

**Output**
--------------


* **individual_images** : Plots each individual trajectory on the first frame of the corresponding film, and saves it as a .jpg file.
* **individual_txt** : Saves the point coordinates of each individual trajectory into a .txt file.
* **group_image** : Plots all trajectories found on a film on its first frame, and saves it as a .jpg file.

|

**Output, transport parameters and statistical analysis**
=============================================================

**Output**
--------------

The main output of the feature finding phase consists of two .csv files : 


* ``filename``.csv : Raw trackpy output containing coordinates of each trajectory.
* ``filename`` _rejoined.csv : Rejoined and filtered trajectories.

Optionally, the script will also generate : 


* A plot of each individual trajectory.
* A .txt files containing the coordinates of each individual trajectory.
* Plots of all trajectories found per film.

The data extraction phase will also generate two .csv files : 


* Per phase parameters.csv : Transport parameters calculated for each phase of each trajectory.
* Trajectory average parameters.csv : Transport parameters averaged from phases of each trajectory.

The statistical analysis phase will generate several files : 


* Barplots for each transport parameters, as .jpg files.
* Boxplots for each transport parameters, as .jpg files.
* Dunn's test tables for each transport parameters, as .jpg files.
* A single .txt file with the p-values for each transport parameters.

|

**Transport parameters**
----------------------------

The following transport parameters are extracted and analyzed from each trajectory.


* **Pausing time** : Time, in seconds, that the feature spent in STOP phases.
* **Pausing frequency** : Frequency at which the feature paused, in number of events per minute.
* **Curvilign velocity** : Also known as segmental velocity, the speed of the feature in µm/s.
* **Processivity** : Time, in seconds, that the feature spent in GO phases.
* **Run length** : Length, in µm, travelled during GO phases.
* **Diagonal size** : Overall length of the trajectory.
* **Fraction of time paused** : Fraction of the time that the feature spent paused.
* 
  **Fraction of moving particles** : Ratio of moving particles to non-moving particles. 

    It is estimated by dividing the number of trajectories analyzed for each file by the number of features found on the first frame of that file. 
    It does not take into account trajectories that were filtered out before analysis, or features that might appear after the first frame. 
    It is therefore not an absolute measure of the fraction of moving particles, and should only be used for relative comparison between experimental conditions.

If the antero_retro setting is enabled : 


* Some of the parameters will be duplicated for anterograde and retrograde transport.
* **Directionality** : ratio of anterograde to retrograde transport. 1 means a purely anterograde transport, 0 a purely retrograde transport.

Additionally : 


* Intensity : Average integrated brightness of the feature over the course of the trajectory. Separated between GO and STOP phases.
* Variance : Standard deviation of the intensity. Similarly separated between GO and STOP phases.
* Number of stops : Total number of pauses within a trajectory.
* Phase-specific parameters :

  * Phase code : 2 signifies a GO phase, 0 a STOP phase.
  * Phase length : Length, in points, of the phase.
  * Vectorial velocity : Speed calculated from the Euclidean distance between the first and last point of the phase.
  * Phase duration : Duration of the phase, in seconds.

|

.. _Statistical analysis:

**Statistical analysis**
----------------------------
This part of the script statistically compares transport parameters for each conditions.

The script first checks for normality of distribution for each parameter. It then applies appropriate statistical tests : 

* First, a Kruskal-Wallis test is applied to check for statistically significant differences between each conditions.
* Then, a post-hoc Dunn's test is applied to check for pair-wise differences.

Barplots and boxplots are generated for each parameter as well.

Results from the Krusal-Wallis as well as normality tests are stored in a single .txt file.

Plots are saved as .jpg files. Dunn's test results are also stored as tables in .jpg files.

A few caveats : 

* Kruskal-Wallis is used even in cases where there are only two conditions. This is because ``scipy.stats``'s Mann-Whitney U test lacks a ``nan_policy``, which interferes with calculations on parameters where some trajectories might lack data (e.g. retrograde transport parameters in a purely anterograde trajectory).
* t-tests and other appropriate parametric tests have yet to be implemented.
* Conditions to be compared are, for now, simply determined by folder structure. 
  
  * For unidirectional transport, folder structure is as such : ``input_folder/experiment/condition/files.tif``
  * For anterograde and retrograde transport, folder structure is a such : ``input_folder/experiment/condition/animal/eye/files.tif``