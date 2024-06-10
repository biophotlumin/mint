**MINT**
========

[<img src="https://zenodo.org/badge/DOI/10.5281/zenodo.5713847.svg">](<(https://doi.org/10.5281/zenodo.5713847)>)

**M**odified **I**ntraneuronal **N**anoparticle **T**racking (**MINT**) is a Python script used to extract intraneuronal transport parameters from trajectories of optically active nanoparticles.

It relies on [**trackpy**](https://github.com/soft-matter/trackpy) to extract trajectories from videos, from which transport parameters are then calculated and statistically tested.

It can reasonably be applied to videos of fluorophore-labelled organelles or similar experiments.

Please refer to the [**documentation**](https://lumin-mint.readthedocs.io/en/latest/) for further information.

Quickstart
---

Install the latest version : 

    git clone https://github.com/biophotlumin/mint.git

    cd mint

    pip install .

You can then run MINT with :

``mint -f <input_folder> -p <path/to/file> -l -e -s``

``-f, --folder`` : Path to raw data folder.

``-p, --params`` : Path to config file (YAML or JSON).

``-l, --locate`` : Run tracking.

``-e, --extract`` : Extract transport parameters from trajectories.

``-s, --stats`` : Run statistical analysis.
