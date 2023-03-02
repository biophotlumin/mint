Install & run
======================

If you are starting from scratch, we recommend installing `Anaconda <https://docs.anaconda.com/anaconda/install/>`_.

Then, download the latest release and extract it into a dedicated folder.

You can create a dedicated environment with ``conda create --name mint python=3.10``, 
and activate it with ``conda activate mint``.

You can then install dependencies with ``pip install -r requirements.txt``.

If you have trouble installing ``cvxpy`` with ``pip``, please use ``conda``.

If ``requirements.txt`` does not install dependencies properly, you can also install them individually with ``pip`` or ``conda``.

You can check installed packages and their versions in the current environment with ``conda list``.

**CLI**
^^^^^^^^^^^
You can run MINT with :

``python mint.py -f <input_folder> -p <path/to/file> -l -e -s``

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

If ``--folder`` isn't specified, the script defaults to the path listed in ``mint.py``. 
If that path is empty too, it then defaults to the current working directory.

If ``--params`` isn't specified, the script defaults to the dictionaries listed in ``mint.py``.

When running ``--extract`` without ``--locate`` first, the input folder must contain the results of a previous ``--locate`` run.

Similarly, when running ``--stats`` without ``--extract`` first, the input folder must contain the results of a previous ``--extract`` run.

If neither ``--locate``, ``--extract`` or ``--stats`` are used, the script will go through a full run and generate an experiment report.

Alternatively, you can directly edit and run ``mint.py``.

|
|

For the sake of reproducibility, you can download the ``paper`` branch from the repo, then 
`import <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file>`_ 
the environment we used in Grimaud *et al.* 2022 from ``reproducible_env.yml``.