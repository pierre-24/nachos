"""
Analyze the results stored in data file
"""

import argparse
import os

from qcip_tools import derivatives, derivatives_e
from qcip_tools.chemistry_files import chemistry_datafile

from nachos import exit_failure
from nachos.core import analyzing, shaking

__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'


def treat_only_arg(only_arg):
    only = []
    for x in only_arg.split(';'):

        try:
            d = derivatives.Derivative(x)
        except derivatives.RepresentationError:
            raise ValueError('wrong derivative: {}'.format(x))

        if derivatives.is_geometrical(d):
            raise ValueError('no vibrational contribution is available for {}'.format(x))

        if d.order() < 1:
            raise ValueError('no vibrational contribution for energy')

        only.append(d)

    return only


to_order = {'m': 1, 'a': 2, 'b': 3, 'g': 4}
coordinates_translate = dict((b, a) for a, b in derivatives.COORDINATES.items())
order_to_class = {
    1: derivatives_e.ElectricDipole,
    2: derivatives_e.PolarisabilityTensor,
    3: derivatives_e.FirstHyperpolarisabilityTensor,
    4: derivatives_e.SecondHyperpolarizabilityTensor
}

need_dipole = ['beta_parallel', 'beta_perpendicular', 'beta_kerr']


def treat_properties_arg(properties_arg, dipole=None):
    properties = {}
    for x in properties_arg.split(';'):
        if not x:
            continue

        info = x.split(':')

        if info[0] not in to_order:
            raise ValueError('{} is not allowed'.format(info[0]))

        order = to_order[info[0]]

        if order not in properties:
            properties[order] = []

        if info[1] == '':
            if len(info) != 3:
                raise ValueError('wrong definition {}'.format(x))
            try:
                component = [coordinates_translate[a] for a in info[2]]
            except KeyError:
                raise ValueError('wrong coordinate definition {}'.format(info[2]))

            if len(component) != order:
                raise ValueError('wrong number of coordinates ({}) for {}'.format(len(component), info[0]))

            properties[order].append(
                analyzing.GetPropertyOfTensor(
                    analyzing.component_access, explain='component ' + info[2], component=tuple(component)))
        else:
            if len(info) != 2:
                raise ValueError('wrong definition {}'.format(x))

            k = order_to_class[order]
            if not hasattr(k, info[1]):
                raise ValueError('object {} has not function to compute {}'.format(info[0], info[1]))

            kw = {
                'accessor': analyzing.property_access,
                'converter': analyzing.get_tensor_converter(k),
                'explain': info[1],
                'function': info[1]
            }

            if order == 3 and info[1] in need_dipole:
                if dipole is None:
                    raise ValueError('requested {}, but no dipole available!'.format(info[1]))
                kw['dipole'] = dipole.components

            properties[order].append(analyzing.GetPropertyOfTensor(**kw))

    return properties


# program options
def get_arguments_parser():
    arguments_parser = argparse.ArgumentParser(description=__doc__)
    arguments_parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    arguments_parser.add_argument(
        '-d', '--data', type=str, help='Input h5 file', default='molecule_nd.h5')

    arguments_parser.add_argument(
        '-O', '--only', help='only fetch the properties of given derivatives')

    arguments_parser.add_argument(
        '-p', '--properties', help='which properties to carry out', required=True)

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
        return exit_failure('error while opening data file: {}'.format(str(e)))

    mu = None
    if 'F' in df.derivatives:
        mu = df.derivatives['F']['static']
        print('(! found a dipole moment)')

    try:
        properties = treat_properties_arg(args.properties, mu)
    except ValueError as e:
        return exit_failure('error while treating properties: {}'.format(str(e)))

    only = None
    if args.only:
        try:
            only = treat_only_arg(args.only)
        except ValueError as e:
            return exit_failure('error while treating derivatives: {}'.format(str(e)))

    if not properties:
        return exit_failure('no property to analyze')

    vibrational_contributions = shaking.load_vibrational_contributions(args.data, df.spacial_dof)
    analyzer = analyzing.Analyzer(df, vibrational_contributions)

    analyzer.analyze(properties, only=only)

if __name__ == '__main__':
    main()
