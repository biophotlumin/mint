********************
**User guide**
********************
.. contents:: :backlinks: None
 

**General information**
===============================

Input data must be in the form of .tif or Nikon .nd2 files, and should be placed in a dedicated folder. Proper folder structure is detailed in `Statistical analysis`_.

Default parameters are provided in ``mint.py``, but as they are very specific to our experimental setup please expect some trial and error
before obtaining satisfactory results. 

You can use the ``test_locate.py`` module to try feature findind parameters.


**Parameters and settings**
===============================


* **extension_in** : File format of the videos. Currently, only .tif and .nd2 are supported.


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


  * **search_range** : Maximum distance a feature can move between frames.
  * **memory** (optional) : Maximum number of frames a feature can disappear for and still be considered the same spot.
  * **adaptive_stop** (optional) : Bottom line ``search_range`` that trackpy will gradually move towards if the original parameter results in an unsolvable calculation.
  * **adaptive_step** (optional) : Factor by which ``search_range`` will decrease until it reaches ``adaptive_stop``.
  * **stub_filtering** (setting, optional) : Filters trajectories based on length.

    * **stub_filtering** (parameter) : Minimum number of points required for a trajectory to be retained.

* 
  **Additionally** : 


  * **MSD** (setting, optional) : Computes the Mean Square Displacement (MSD) for each trajectory, and filters them accordingly.

    * **msd** (parameter) : MSD value below which a trajectory is discarded.

  * **rejoining** (optional) : Rejoins trajectories that weren't linked by trackpy in the first place. This function is kept deliberately stringent (i.e. each trajectory cannot be rejoined more than once) to avoid aberrant trajectories and false positives.

    * **threshold_t** : Maximum number of frames between the first and last points of two trajectories for them to be considered for rejoining.
    * **threshold_r** : Maximum distance (in pixels) between the first and last points of two trajectories for them to be considered for rejoining.

  * **SNR_estimation** (optional): Calcualtes the signal-to-noise ratio.

    * **base_level** : Base noise level.

|

**Data extraction**
-----------------------

This part of the script calculates transport parameters from extracted trajectories.


* **r_conf_cut** : Cutoff level for confinement ratio above which a point is considered in a GO phase.
* **px** : Pixel size in µm.
* **dt** : Sampling period, in seconds, as in the amount of time between two frames.
* **min_thr_prec** : Minimum theoretical precision, in nm.
* **sliding_window** : Size, in frames, of the sliding window along which the confinement ratio is calculated.
* **polynomial_fit** (optional) : Filters trajectories based on how well they fit to a 3rd degree polynom.

  * **threshold_poly3** : Tolerated deviation from the 3rd degree polynom.
  * **len_cutoff** : Size, in points, below which trajectories are eliminated.

* **minimization** (optional) : Experimental function of convex minimization of acceleration used to smooth trajectories with a lot of spatial noise.

  * **sigma** : Estimated noise, as the standard deviation of the precision of localization, in nm.

* **antero_retro** (optional) : Separates transport parameters into anterograde and retrograde categories. Note : this is highly dependent on our original experimental setup and will most likely not work elsewhere. 

* **conf_list** (optional) : Saves the confinement ratio of each point into a .csv file that contain a list of points for each trajectory. Can be used to determinate ``r_conf_cut``.

|

**Output**
--------------


* **individual_images** (optional) : Plots each individual trajectory on the first frame of the corresponding film, and saves it.
* **individual_txt** (optional) : Saves the point coordinates of each individual trajectory into a .txt file.
* **group_image** (optional) : Plots all trajectories found on a film on its first frame, and saves it.
* **ordering** (optional) : Specify the order of experimental conditions in graphs.

  * **order** : List of experimental conditions.
* **extension_out** : File format under which graphs will be saved. Can be anything ``matplotlib`` supports.
* **dpi** (optional if ``extension_out`` is vectorial) : DPI of the saved graphs for non-vectorial file formats.
* **clean_up** (optional) : Wether or not to delete individual graph files once they've been included in the experiment report.

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


* ``Per phase parameters.csv`` : Transport parameters calculated for each phase of each trajectory.
* ``Trajectory average parameters.csv`` : Transport parameters averaged from phases of each trajectory.

The statistical analysis phase will generate several files : 


* **Barplots** for each transport parameters.
* **Boxplots** for each transport parameters.
* A single .txt file with the p-values for each transport parameters as well as some other statistics.

Additionally, several dictionaries are dumped as .txt files : 


* ``log.txt`` contains some information about the run.
* ``parameters.txt`` lists the parameters that were used.
* ``settings.txt`` lists the settings that were used.
* ``var.txt`` lists the variables statistically tested.

In the case of a full run, the script will also generate a complete experiment report into a .pdf file.

|

**Transport parameters**
----------------------------

The following transport parameters are extracted and analyzed from each trajectory.


* **Pausing time** : Time, in seconds, that the particle spent in STOP phases.
* **Pausing frequency** : Frequency at which the particle paused, in number of events per minute.
* **Curvilign velocity** : Also known as segmental velocity, the speed of the particle in µm/s.
* **Processivity** : Time, in seconds, that the particle spent in GO phases.
* **Run length** : Length, in µm, travelled during GO phases.
* **Diagonal size** : Overall length of the trajectory.
* **Fraction of time paused** : Fraction of the time that the particle spent paused.
* 
  **Fraction of moving particles** : Ratio of moving particles to non-moving particles. 

    It is estimated by diving the number of trajectories analyzed for each file by the number of features found on the first frame of a film. 

    It does not take into account trajectories that were filtered out before analysis, or features that might appear after the first frame. 

    It is therefore not an absolute measure of the fraction of moving particles, and should only be used for relative comparison between experimental conditions.

If the antero_retro setting is enabled : 


* Some of the parameters will be duplicated for anterograde and retrograde transport.
* **Directionality** : ratio of retrograde to anterograde transport. 1 means a purely retrograde transport, 0 a purely anterograde transport.

Additionally : 


* **Intensity** : Average integrated brightness of the feature over the course of the trajectory. Separated between GO and STOP phases.
* **Variance** : Standard deviation of the intensity. Similarly separated between GO and STOP phases.
* **Number of stops** : Total number of pauses within a trajectory.
* Phase-specific parameters :

  * **Phase code** : 2 signifies a GO phase, 0 a STOP phase.
  * **Phase length** : Length, in points, of the phase.
  * **Vectorial velocity** : Speed calculated from the Euclidean distance between the first and last point of the phase.
  * **Phase duration** : Duration of the phase, in seconds.

|

**Statistical analysis**
----------------------------
This part of the script statistically compares transport parameters between each experimental condition.

The script first checks for normality of distribution for each parameter. It then applies appropriate statistical tests : 

* If there are two experimental conditions and the distribution is normal, a Student's t-test is applied. If it is not normal, a ranksums test is applied.
* If there are more than two experimental conditions, a Kruskal-Wallis test is applied. Then, a post-hoc Dunn's test is applied to check for pair-wise differences.

Barplots and boxplots are generated for each parameter as well.

Results from the statistical tests are stored in a single .txt file.

* Conditions to be compared are, for now, simply determined by folder structure, such as :

 ``input_folder/experiment/condition/replicate/sample/files.tif``