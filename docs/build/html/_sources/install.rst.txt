Install
======================

If you are starting from scratch, we recommend installing `Anaconda <https://docs.anaconda.com/anaconda/install/>`_.

Then, download the latest release and extract it into a dedicated folder.

You can then create an environment with proper dependencies by using :

``conda create --name <env_name> --file requirements.txt``, 


or install dependencies to an already existing environment with :

``pip install -r requirements.txt``.

If you have trouble installing ``cvxpy`` with ``pip``, please use ``conda``.

If ``requirements.txt`` does not install dependencies properly, you can install them individually with ``pip`` or ``conda``.

**CLI**
^^^^^^^^^^^
Run ``cli.py``.

**GUI**
^^^^^^^^^^^
Run ``gui.py``.

**Script**
^^^^^^^^^^^
Edit the ``parameters`` and ``settings`` dictionaries found in ``script.py``, then run it either through command line or directly from an IDE.

|
|

For the sake of reproducibility, you can otherwise download the ``publication`` release, then 
`import <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file>`_ 
the environment we used for ``publication`` from ``environment.yml``.