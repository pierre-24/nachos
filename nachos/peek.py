#!/usr/bin/env python3
"""
Peek into an output file (LOG, FCHK, H5, ...)
"""

import argparse
import sys
from scipy import constants

from qcip_tools import derivatives_e, quantities, derivatives_g
from qcip_tools.chemistry_files import helpers, PropertyNotPresent

import nachos

__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'


def to_nanometer(val):
    """Convert frequency to nanometer

    :param val:
    :return:
    """

    if type(val) is str:
        return val

    converted = \
        constants.h * constants.c / (val * quantities.convert(quantities.ureg.hartree, quantities.ureg.joule)) * 1e9

    return '{:.1f}nm'.format(converted)


def print_electric_d_tensor(electrical_derivatives, representation):
    if representation in electrical_derivatives:
        freqs = [x for x in electrical_derivatives[representation].keys()]
        freqs.sort(key=lambda x: derivatives_e.convert_frequency_from_string(x))

        name = derivatives_e.NAMES[representation]

        # include dipole if found
        kw = {}
        if len(representation) == 3:
            if 'F' in electrical_derivatives:
                kw['dipole'] = electrical_derivatives['F']['static'].components

        for freq in freqs:
            print('{}, w={} ({:.6f} a.u.)'.format(
                name, to_nanometer(freq), derivatives_e.convert_frequency_from_string(freq)))

            print(electrical_derivatives[representation][freq].to_string(**kw))


def print_geometric_d_tensor(geometrical_derivatives, representation, molecule):
    if representation in geometrical_derivatives:
        name = derivatives_g.NAMES[representation]
        print(name)
        if molecule is not None:
            print(geometrical_derivatives[representation].to_string(molecule=molecule))
        else:
            print(geometrical_derivatives[representation].to_string())


def get_arguments_parser():
    arguments_parser = argparse.ArgumentParser(description=__doc__)
    arguments_parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    arguments_parser.add_argument(
        'infile',
        nargs='?',
        default=sys.stdin,
        action=helpers.create_open_chemistry_file_action(),
        help='source of the derivatives')

    return arguments_parser


def main():
    args = get_arguments_parser().parse_args()

    print('This is nachos_peek (v{}) from nachos (v{})\n'.format(__version__, nachos.__version__))

    to_show = [
        # alpha
        'FF', 'dD',
        # beta
        'FFF', 'dDF', 'FDd', 'XDD',
        # gamma
        'FFFF', 'dDFF', 'dFFD', 'XDDF', 'dDDd', 'XDDD',

        # cartesian
        'G', 'GG', 'GGG'
        # normal
        'N', 'NN', 'NNN'
    ]

    # (pure) electrical derivatives
    if args.infile.has_property('electrical_derivatives'):
        try:
            electrical_derivatives = args.infile.property('electrical_derivatives')

            if 'F' in electrical_derivatives:
                print('dipole moment:')
                print(electrical_derivatives['F']['static'].to_string())
                print('')

            for f in to_show:
                print_electric_d_tensor(electrical_derivatives, f)
        except PropertyNotPresent:
            pass

    #  (pure) geometrical derivatives
    if args.infile.has_property('geometrical_derivatives'):

        try:
            molecule = args.infile.property('molecule')
        except PropertyNotPresent:
            molecule = None

        try:
            geometrical_derivatives = args.infile.property('geometrical_derivatives')

            for f in to_show:
                print_geometric_d_tensor(geometrical_derivatives, f, molecule)
        except PropertyNotPresent:
            pass


if __name__ == '__main__':
    main()
