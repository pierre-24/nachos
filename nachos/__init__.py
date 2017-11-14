"""
NACHOS: numerical differentiation code
"""
import argparse
import os

__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'

# List of the scripts that are installed, without the .py extension. The name is used to give the command name.
provide_scripts = [
    'nachos_prepare',
    'nachos_cook'
]


def make_console_scripts(package_name='nachos'):
    """Function used to generate the ``console_scripts`` list for ``setup.py``

    :rtype: list
    """

    console_scripts = []
    for script in provide_scripts:
        path = os.path.join(package_name, script + '.py')
        if not os.path.isfile(path):
            raise FileNotFoundError(path)

        console_scripts.append('{0} = {1}.{0}:main'.format(script, package_name))

    return console_scripts


def is_dir(dirname):
    """Checks if a path is an actual directory"""
    if not os.path.isdir(dirname):
        msg = '{0} is not a directory'.format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname
