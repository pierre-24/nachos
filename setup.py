import nachos

from setuptools import setup
try: # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError: # for pip <= 9.0.3
    from pip.req import parse_requirements

pkgs = []
dependency_links = []
for pkg in parse_requirements('requirements.txt', session=False):
    if pkg.req:
        pkgs.append(str(pkg.req))

setup(
    name='nachos',
    packages=['nachos', 'nachos.core'],
    version=nachos.__version__,
    author=nachos.__author__,
    author_email=nachos.__email__,
    description=nachos.__doc__,
    classifiers=[
        'Environment :: Scientific',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3'
    ],
    install_requires=pkgs,
    dependency_links=dependency_links,
    python_requires='>=3',
    test_suite='tests',
    entry_points={
        'console_scripts': nachos.make_console_scripts()
    },
)
