from setuptools import setup, find_packages
from os import path

import nachos

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md')) as f:
    long_description = f.read()

with open(path.join(here, 'requirements/requirements-base.in')) as f:
    requirements = f.readlines()

with open(path.join(here, 'requirements/requirements.in')) as f:
    requirements_dev = f.readlines()[1:]

setup(
    name='nachos',
    version=nachos.__version__,

    # Description
    description=nachos.__doc__,
    long_description=long_description,
    long_description_content_type='text/markdown',
    keywords='website',

    project_urls={
        'Bug Reports': 'https://github.com/pierre-24/nachos/issues',
        'Source': 'https://github.com/pierre-24/nachos',
    },

    url='https://github.com/pierre-24/nachos',
    author=nachos.__author__,

    # Classifiers
    classifiers=[
        'Environment :: Scientific',
        'Operating System :: OS Independent',

        # Specify the Python versions:
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9'
    ],

    packages=find_packages(),
    python_requires='>=3.8',
    test_suite='tests',
    entry_points={
        'console_scripts': nachos.make_console_scripts()
    },

    # requirements
    install_requires=requirements,

    extras_require={  # Optional
        'dev': requirements_dev,
    },
)