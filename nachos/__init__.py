"""
NACHOS: numerical differentiation code
"""
import argparse
import os
import sys

__version__ = '0.2.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'

# List of the scripts that are installed, without the .py extension. The name is used to give the command name.
provide_scripts = [
    'make',
    'prepare',
    'cook',
    'bake',
    'shake',
    'analyze'
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

        console_scripts.append('{0} = {1}.{2}:main'.format('nachos_' + script, package_name, script))

    return console_scripts


def is_dir(dirname):
    """Checks if a path is an actual directory"""
    if not os.path.isdir(dirname):
        msg = '{0} is not a directory'.format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname


def exit_failure(msg, status=1):
    """Write a message in stderr and exits

    :param msg: the msg
    :type msg: str
    :param status: exit status (!=0)
    :type status: int
    """

    sys.stderr.write(msg)
    sys.stderr.write('\n')
    return sys.exit(status)
