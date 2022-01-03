from setuptools import setup, find_packages

setup(
    name='cropsr',
    version='0.1.0',
    description='CROPSR package for genome-wide CRISPR gRNA analysis and evaluation',
    long_description=open('README.md').read(),
    url='https://github.com/H2muller/CROPSR/',
    author='Hans Müller Paul',
    author_email='hmpaul2@illinois.edu'
)

"""
License=APACHE-2.0
classifiers=[
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'License :: OSI Approved :: Apache Software License',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9'
]
"""

keywords='CROPSR genome-wide CRISPR gRNA analysis evaluation crop genomics knockout'
packages=find_packages(exclude=['docs', 'tests*'])
install_requires=['pandas argparse']
package_data={
    'sample': ['package_data.dat'],
}
data_files=None
