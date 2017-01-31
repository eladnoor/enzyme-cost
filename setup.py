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
    install_requires=['scipy>=0.18.1',
                      'numpy>=1.12.0',
                      'pulp>=1.6.1',
                      'pyparsing>=2.1.10'],
    data_files=[('ecm_ecoli_aerobic', ['data/ecm_ecoli_aerobic.tsv']),
               ],
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2.7',
    ],
)

