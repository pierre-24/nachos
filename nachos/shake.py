"""
Shake it! (compute the vibrational contributions)
"""

import argparse
import os

from qcip_tools import derivatives
from qcip_tools.chemistry_files import chemistry_datafile

from nachos import exit_failure
from nachos.core import shaking

__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'


def treat_only_arg(only_arg):
    only = []
    for x in only_arg.split(';'):
        info = x.split(':')
        level = 2

        if len(info) > 2:
            raise ValueError('{} is not a correct input'.format(x))

        try:
            d = derivatives.Derivative(info[0])
        except derivatives.RepresentationError:
            raise ValueError('wrong derivative: {}'.format(info[0]))

        if derivatives.is_geometrical(d):
            raise ValueError('no vibrational contribution is available for {}'.format(info[0]))

        if d.order() < 1:
            raise ValueError('no vibrational contribution for energy')

        if len(info) > 1:
            try:
                level = int(info[1])
            except ValueError:
                raise ValueError('wrong level {} for {}'.format(info[1], info[0]))

        if level < 0:
            raise ValueError('wrong level {} for {} (should be positive)'.format(level, info[0]))

        only.append((d, level))

    return only


def treat_frequencies_arg(frequencies_arg):
    frequencies = []
    for x in frequencies_arg.split(';'):
        allowed_units = ['nm', 'ev', 'cm-1']

        value = None

        try:
            value = float(x)
        except ValueError:
            for unit in allowed_units:
                if x[-len(unit):] == unit:
                    try:
                        float(x[:-len(unit)])
                        value = x
                    except ValueError:
                        pass
                    break

        if value is None:
            raise ValueError('{} is not an allowed frequency'.format(x))

        frequencies.append(value)

    return frequencies


# program options
def get_arguments_parser():
    arguments_parser = argparse.ArgumentParser(description=__doc__)
    arguments_parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    arguments_parser.add_argument(
        '-d', '--data', type=str, help='Input h5 file', default='molecule_nd.h5')

    arguments_parser.add_argument('-V', '--verbose', type=int, help='Level of details (0 to 3)', default=0)

    arguments_parser.add_argument(
        '-O', '--only', help='only compute the contribution of given derivatives')

    arguments_parser.add_argument(
        '-F', '--frequencies', help='compute vibrational contribution for set of frequencies')

    return arguments_parser


# main
def main():
    args = get_arguments_parser().parse_args()

    df = chemistry_datafile.ChemistryDataFile()

    if not os.path.exists(args.data):
        return exit_failure('data file {} does not exists'.format(args.data))

    try:
        with open(args.data) as f:
            df.read(f)
    except chemistry_datafile.BadChemistryDataFile as e:
        return exit_failure('error whil opening data file: {}'.format(str(e)))

    shaker = shaking.Shaker(df)

    only = None
    if args.only:
        try:
            only = treat_only_arg(args.only)
        except ValueError as e:
            return exit_failure('error while treating arg: {}'.format(str(e)))

    frequencies = None
    if args.frequencies:
        try:
            frequencies = treat_frequencies_arg(args.frequencies)
        except ValueError as e:
            return exit_failure('error while treating frequencies: {}'.format(str(e)))

    try:
        shaker.shake(verbosity_level=args.verbose, only=only, frequencies=frequencies)
    except shaking.BadShaking as e:
        return exit_failure('error while making recipe: {}'.format(str(e)))

if __name__ == '__main__':
    main()
