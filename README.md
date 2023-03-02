**MINT**
========

**M**odified **I**ntraneuronal **N**anoparticle **T**racking (**MINT**) is a Python script used to extract intraneuronal transport parameters from the trajectories of optically active nanoparticles.

It relies on [**trackpy**](https://github.com/soft-matter/trackpy) to extract trajectories from video files, from which transport parameters are then calculated and statistically tested.

Its purpose is to automate analysis to the point where an input of raw video files results in an output of graphs and p-values.

It can reasonably be applied to videos of fluorophore-labelled organelles, or any object that can be tracked with trackpy.

Please refer to the [**documentation**](https://lumin-mint.readthedocs.io/en/latest/) for further information.

Quickstart
---

Download the lastest release and extract its contents to a dedicated folder, then install the dependencies found in `requirements.txt`.

You can then run MINT with :

``python mint.py -f <input_folder> -p <path/to/file> -l -e -s``

``-f, --folder``
    Path to raw data folder.

``-p, --params``
    Path to config file (YAML or JSON).

``-l, --locate``
    Run tracking.

``-e, --extract``
    Extract transport parameters from trajectories.

``-s, --stats``
    Run statistical analysis.
