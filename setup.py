import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = 'Ensembler2',
    version = '0.2',
    author = 'Patrick B. Grinaway',
    author_email = 'patrick.grinaway@choderalab.org',
    description = 'Generation of diverse protein structural ensembles, for the initialization of molecular dynamics simulations and subsequent construction of Markov state models. ',
    long_description = read('README.md'),
    packages = ['Ensembler2', 'tests'],
    data_files = [('', ['LICENSE']), ('templates', ['project-data.yaml-TEMPLATE', 'manual-specifications.yaml-TEMPLATE'])],
)
