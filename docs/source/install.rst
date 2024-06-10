Install & run
======================

If you are starting from scratch, we recommend installing `Anaconda <https://docs.anaconda.com/anaconda/install/>`_.

MINT has been tested with Python 3.10, 3.11 and 3.12.

Then, download the latest release and extract it into a dedicated folder, or clone the repo with : 

.. code:: console

    $ git clone https://github.com/biophotlumin/mint.git
    $ cd mint

Then create a dedicated environment with : 

.. code:: console

    $ conda create --name mint python=3.10
    $ conda activate mint

Then install MINT with :

.. code:: console

    $ pip install .

You can also install optional dependencies by passing one of the following as : 

.. code:: console

    $ pip install .[optional]

``notebook``
    Jupyter dependencies to use interactive notebooks.

``bioformats``
    ImageJ interface to read Bio-Formats extensions.

``gui``
    GUI dependencies.

``solvers``
    CVXPY and associated solvers for experimental trajectory denoising.

``all``
    All of the above.

If you have trouble installing ``cvxpy`` with ``pip``, try using ``conda``.

For development, you can install in editable mode :

.. code:: console

    $ pip install -e .[optional]

You can then install :

``docs``
    Sphinx and other dependencies to build docs.

``dev``
    All of the above, including ``docs``.

**CLI**
^^^^^^^^^^^
You can run MINT with :

.. code:: console

    $ mint -f <input_folder> -p <path/to/file> -l -e -s

``-f, --folder``
    Optional. Path to raw data folder.

``-p, --params``
    Optional. Path to config file (YAML or JSON).

``-l, --locate``
    Optional. Run tracking.

``-e, --extract``
    Optional. Extract transport parameters from trajectories.

``-s, --stats``
    Optional. Run statistical analysis.

Ideally, the input folder should be specified in the config file (see the `user guide <user_guide.html>`_).

If the config file does not contain an input path and ``--folder`` isn't specified, the script will default to the current working directory.

If only ``--folder`` is specified, the script will look for a config file in that folder. Otherwise, it will fall back to default parameters.

When running ``--extract`` without ``--locate`` first, the input folder must contain the results of a previous ``--locate`` run.

Similarly, when running ``--stats`` without ``--extract`` first, the input folder must contain the results of a previous ``--extract`` run.

If neither ``--locate``, ``--extract`` or ``--stats`` are specified, the script will go through a full run.

|

For the sake of reproducibility, you can download the ``paper`` branch from the repo, then 
`import <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file>`_ 
the environment we used in Grimaud *et al.* 2022 from ``reproducible_env.yml``.