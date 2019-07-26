from setuptools import setup, find_packages

with open('README.rst', mode='r') as f:
    l_description = f.read()

setup(name='randmag',
      version='1.2.0',
      description='Construct artifical MAGs from complete genomes according to distribution',
      long_description=l_description,
      url='https://bitbucket.org/EricHugo/randmag',
      author='Eric Hugoson',
      author_email='eric@hugoson.org',
      license='GNU General Public License v3 or later (G  PLv3+)',
        
      classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Environment :: Console',
            'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
        ],
      keywords='bioinformatics simulation metagenomics',
      install_requires=[
            'biopython>=1.72',
            'numpy>=1.16.1',
            'matplotlib>=2.0.2',
            'pandas>=0.23.4',
            'seaborn>=0.9.0',
            'scipy>=1.1.0',
            ],
      packages=find_packages(),
      include_package_data=True,
      entry_points={
            'console_scripts': [
                'randMAG = randmag.randmag:main',
                ],
            }
      )

