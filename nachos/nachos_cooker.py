"""
Create the input files for the quantum chemistry programs
"""

import argparse

__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'


# program options
def get_arguments_parser():
    arguments_parser = argparse.ArgumentParser(description=__doc__)
    arguments_parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    return arguments_parser


# main
def main():
    args = get_arguments_parser().parse_args()


if __name__ == '__main__':
    main()
