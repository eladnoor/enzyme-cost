from setuptools import setup

setup(
    name = 'EnzymeCost',
    version = '0.1.0',
    author = 'Elad Noor',
    author_email='noor@imsb.biol.ethz.ch',
    description = 'A collection of tools that perform Enzyme Cost Minimization',
    license = 'MIT',
    packages=['ecm'],
    url='https://github.com/eladnoor/enzyme-cost',
    install_requires=['scipy>=1.0.0',
                      'numpy>=1.13.3',
                      'pulp>=1.6.8',
                      'matplotlib>=2.1.0',
                      'pandas>=0.21.0',
                      'tablib>=0.12.1',
                      'equilibrator_api>=0.1.0'],
    data_files=[('ecm_ecoli_aerobic', ['data/ecoli_ccm_aerobic_ProteinComposition_haverkorn_ECM_Model.tsv']),
               ],
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3.6',
    ],
)

