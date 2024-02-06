from setuptools import setup

setup(
    name='mint',
    version='0.3.0',    
    description='Modified Intraneuronal Nanoparticle Tracking (MINT) is a Python script used to extract intraneuronal transport parameters from the trajectories of optically active nanoparticles.',
    url='https://github.com/biophotlumin/mint',
    author='Baptiste Grimaud',
    author_email='baptiste.grimaud@ens-paris-saclay.fr',
    license='GNU General Public License v3.0',
    packages=['mint'],
    install_requires=['trackpy>=0.6.0',
                      'numpy<=1.26',
                      'numba',
                                           
                      ],

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Image Processing',
    ],
)