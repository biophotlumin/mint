Install
======================

If you are starting from scratch, we recommend installing `Anaconda <https://docs.anaconda.com/anaconda/install/>`_.

Then, download the latest release and extract it into a dedicated folder.

You can then create an environment with proper dependencies by using ``pip install -r requirements.txt`` or ``conda create --name <env_name> --file requirements.txt``.

**CLI**
^^^^^^^^^^^
Edit the ``parameters`` and ``settings`` dictionaries found in ``script.py``, then run it either through command line or directly from an IDE.

**GUI** (work in progress)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Run ``gui.py``.

|
|

For the sake of reproducibility, you can otherwise download the ``paper`` release, then 
`import <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file>`_ 
the environment we used for ``paper`` from ``environment.yml``.