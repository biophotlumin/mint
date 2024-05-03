********************
User guide
********************

.. contents:: :backlinks: None
 

**Input**
===============================


All parameters and settings are specified in a single ``config.yml`` file.
This file can be moved and renamed.

``input_folder`` is the absolute path to the folder containing videos to analyse.
Proper folder structure is detailed in `Statistical analysis`_.

MINT can read .tif and Nikon .nd2 files through ``imageio`` and ``nd2``, respectively.

Other file types supported by Bio-Formats can be read through ``imagej``, albeit much slower.

Default parameters are provided in ``config.yml`` but are very much experiment-dependent.
Expect some trial and error before obtaining satisfactory results. 

You can use the ``test_locate`` function to optimize feature findind parameters.


**Parameters and settings**
===============================


* **extension_in** : File format of the videos (as its extension, without period).

* **parallel** : ``joblib``-based parallelization. Depending on the number of avaible cores, it can dramatically speed up calculation.

* **parallel_tracking** : Experimental. Parallelizes the entire tracking process. Reading files from disk and available memory can become bottlenecks.


**Preprocessing**
---------------------

Two preprocessing filters can be applied to the videos before going through tracking :


* **Wavelet** denoising reduces background noise and improves the signal-to-noise ratio.


* **Top-hat** transform (white) removes large spots and artifacts.

|

**Feature finding and linking**
-----------------------------------

This part is handled by `\trackpy <https://github.com/soft-matter/trackpy>`_. 
Please refer to their own documentation for more details on each parameter.


* 
  **Feature finding** :


  * **diameter** : Size of the features to be found (in pixels). Must be odd.
  * **minmass** : Minimum integrated brightness.
  * **separation** : Minimum separation between features (in pixels).

* 
  **Linking** : 


  * **search_range** : Maximum distance a feature can move between frames.
  * **memory** (optional) : Maximum number of frames a spot can disappear for and still be considered the same feature.
  * **adaptive_stop** (optional) : Bottom line ``search_range`` that trackpy will gradually move towards if the original parameter results in an unsolvable calculation.
  * **adaptive_step** (optional) : Factor by which ``search_range`` will decrease until it reaches ``adaptive_stop``.
  * **stub_filtering** (setting, optional) : Filters trajectories based on length.

    * **stub_filtering** (parameter) : Minimum number of points required for a trajectory to be retained.

* 
  **Postprocessing** : 


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


* **r_conf_cut** : Cutoff level for confinement ratio above which a point is considered in a GO phase. See the ``rconf`` notebook for more details.
* **px** : Pixel size, in µm.
* **dt** : Sampling period, in seconds.
* **min_thr_prec** : Minimum theoretical precision, in nm.
* **sliding_window** : Size, in frames, of the sliding window along which the confinement ratio is calculated.
* **polynomial_fit** (optional) : Filters trajectories based on how well they fit to a 3rd degree polynom.

  * **threshold_poly3** : Tolerated deviation from the 3rd degree polynom.
  * **len_cutoff** : Size, in points, below which trajectories are eliminated.

* **minimization** (optional) : Experimental function of convex minimization of acceleration used to smooth trajectories with a lot of spatial noise.

  * **sigma** : Estimated noise, as the standard deviation of the precision of localization, in nm.

* **antero_retro** (optional) : Separates transport parameters into anterograde and retrograde categories. Note : this is highly dependent on our original experimental setup and will most likely not work elsewhere. 

* **conf_list** (optional) : Saves the confinement ratio of each point into a .csv file that contain a list of points for each trajectory. Can be used to set ``r_conf_cut``.

* **theta** (optional) : Calculates the orientation of particles based on the variation of their intensity. Only applies to nanoparticles with non-isotropic emission.

|

**Output**
--------------


* **individual_images** (optional) : Plots each individual trajectory on the first frame of the corresponding video, and saves it.
* **individual_txt** (optional) : Saves the point coordinates of each individual trajectory into a .txt file.
* **group_image** (optional) : Plots all trajectories found on a film on its first frame, and saves it.
* **ordering** (optional) : Specify the order of experimental conditions in graphs.

  * **order** : List of experimental conditions.
* **extension_out** : File format under which graphs will be saved. Can be anything ``matplotlib`` supports.
* **dpi** (optional if ``extension_out`` is vectorial) : DPI of the saved graphs for non-vectorial file formats.

|

**Output, transport parameters and statistical analysis**
=============================================================

**Output**
--------------

The main output of the feature finding phase consists of two .csv files per video, placed in folders matching the hierarchy of the input : 


* ``filename``.csv : Raw trackpy output containing coordinates of each trajectory.
* ``filename`` _rejoined.csv : Rejoined and filtered trajectories.

Optionally, the script can also generate : 


* A plot of each individual trajectory.
* A .txt files containing the coordinates of each individual trajectory.
* Plots of all trajectories found per video.

The data extraction phase will also generate two .csv files, placed in a separate folder : 


* ``phase_parameters.csv`` : Transport parameters calculated for each phase of each trajectory.
* ``trajectory_parameters.csv`` : Transport parameters averaged from phases of each trajectory.

Optionally, this folder will also contain a ``Confinement ratio.csv`` file.

The statistical analysis phase will generate several files : 


* **Barplots** for each transport parameters.
* **Boxplots** for each transport parameters.
* A single .txt file with the p-values for each transport parameters as well as some other statistics.

Additionally, several dictionaries are dumped as .yml files : 


* ``log.yml`` contains some information about the run.
* ``parameters.yml`` lists the parameters that were used.
* ``settings.yml`` lists the settings that were used.
* ``vars.yml`` lists the variables statistically tested.

|

**Transport parameters**
----------------------------

The following transport parameters are extracted and analyzed from each trajectory.


* **Pausing time** : Time, in seconds, that the particle spent in STOP phases.
* **Pausing frequency** : Frequency of STOP phases, in number of events per minute.
* **Duration** : Duration of the trajectory, in seconds.
* **Curvilign velocity** : Also known as segmental velocity, the speed of the particle in µm/s.
* **Processivity** : Time, in seconds, that the particle spent in GO phases.
* **Run length** : Length, in µm, travelled during GO phases.
* **Diagonal length** : Distance between the first and last points of the trajectory.
* **Curvilign length** : Sum of all run lengths.
* **Fraction of time paused** : Fraction of the time that the particle spent paused.
* 
  **Fraction of moving particles** : Ratio of moving particles to non-moving particles. 

    It is estimated by diving the number of trajectories analyzed for each file by the number of features found on the first frame of a film. 

    It does not take into account trajectories that were filtered out before analysis, or features that might appear after the first frame. 

    It is therefore not an absolute measure of the fraction of moving particles, and should only be used for relative comparison between experimental conditions.

If the antero_retro setting is enabled : 


* Some of the parameters will be duplicated for anterograde and retrograde transport.
* **Directionality** : Ratio of retrograde to anterograde transport. 1 means a purely retrograde transport, 0 a purely anterograde transport.
* **Switch** : Amount of directionality reversals, i.e. the number of STOP phases in between GO phases of opposite directionality.
* **Switch A to R** : Reversals from anterograde to retrograde.
* **Switch R to A** : Reversals from retrograde to anterograde.
* **Normalized switch** : Amount of reversals normalized to the duration of the trajectory.
* **Pausing time switch** : Pausing time between GO phases of opposite directionality.
* **Pausing time antero** : Pausing time between anterograde GO phases.
* **Pausing time retro** : Pausing time between retrograde GO phases.


Additionally : 


* **Intensity** : Average integrated brightness of the feature over the course of the trajectory. Separated between GO and STOP phases.
* **Variance** : Average standard deviation of the intensity. Similarly separated between GO and STOP phases.
* **Number of stops** : Total number of pauses within a trajectory.
* **Theta** : Variation of the angle of the nanoparticle.
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