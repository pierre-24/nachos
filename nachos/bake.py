"""
From h5 file, perform numerical differentiation
"""

import os
import argparse

from qcip_tools import derivatives_g, derivatives
from qcip_tools.chemistry_files import helpers, PropertyNotDefined, PropertyNotPresent

import nachos
from nachos.core import files, baking
from nachos import exit_failure

__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'


def treat_only_arg(recipe, t):
    """Treat the --only argument (delayed because it needs the recipe to create the ``Derivative`` object)

    :param recipe: the recipe
    :type recipe: nachos.core.files.Recipe
    :param t: the parameter
    :type t: str
    :rtype: list
    """

    only = []
    for x in t.split(';'):
        if not x:
            continue

        info = x.split(':')
        level = 0

        if len(info) > 2:
            raise ValueError('{} is not a correct input'.format(x))
        if len(info) == 2:
            b_repr = info[0]
            try:
                level = int(info[1])
            except ValueError:
                raise ValueError(
                    'level of differentiation {} (associated with {}) is not a number'.format(info[1], info[0]))
        else:
            b_repr = info[0]

        try:
            b = derivatives.Derivative(b_repr if b_repr != 'energy' else '', spacial_dof=recipe.dof)
        except derivatives.RepresentationError:
            raise ValueError('derivative {} is incorrect'.format(b_repr))

        only.append((b, level))

    return only


# program options
def get_arguments_parser():
    arguments_parser = argparse.ArgumentParser(description=__doc__)
    arguments_parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    arguments_parser.add_argument(
        '-r', '--recipe', type=argparse.FileType('r'), help='Recipe file', default='./nachos_recipe.yml')
    arguments_parser.add_argument(
        '-d', '--data', type=str, help='H5 data file (output of nachos_cook)', default='nachos_data.h5')
    arguments_parser.add_argument(
        '-o', '--output', type=str, help='Output h5 file', default='molecule_nd.h5')

    arguments_parser.add_argument('-V', '--verbose', type=int, help='Level of details (0 to 3)', default=0)

    arguments_parser.add_argument(
        '-S', '--do-not-steal', action='store_true', help='do not add base derivatives to output file')

    arguments_parser.add_argument(
        '-O', '--only', help='only differentiate a subset of the original recipe')

    arguments_parser.add_argument(
        '-p', '--project', action='store_true', help='project geometrical derivatives over normal mode (if hessian)')

    arguments_parser.add_argument(
        '-H', '--hessian',
        help='consume hessian from an other file',
        action=helpers.create_open_chemistry_file_action())

    return arguments_parser


# main
def main():
    args = get_arguments_parser().parse_args()

    if args.verbose >= 1:
        print('This is nachos_bake (v{}) from nachos (v{})\n'.format(__version__, nachos.__version__))

    recipe_directory = os.path.dirname(args.recipe.name)
    recipe = files.Recipe(directory=recipe_directory)

    try:
        recipe.read(args.recipe)
    except files.BadRecipe as e:
        return exit_failure('error while opening recipe: {}'.format(str(e)))

    if not os.path.exists(os.path.join(recipe_directory, args.data)):
        return exit_failure('Data file {} does not exists'.format(os.path.join(recipe_directory, args.data)))

    storage = files.ComputationalResults(recipe, directory=recipe_directory)
    storage.read(args.data)

    baker = baking.Baker(recipe, storage, directory=recipe_directory)
    only = None

    if args.only:
        try:
            only = treat_only_arg(recipe, args.only)
        except ValueError as e:
            return exit_failure('error while treating --only: {}'.format(str(e)))

    try:
        cf = baker.bake(copy_zero_field_basis=not args.do_not_steal, verbosity_level=args.verbose, only=only)
    except baking.BadBaking as e:
        return exit_failure('error while baking: {}'.format(str(e)))

    if args.project:
        hessian = None
        if args.hessian:
            try:
                geometrical_derivatives = args.hessian.property('geometrical_derivatives')
                if 'GG' in geometrical_derivatives:
                    hessian = geometrical_derivatives['GG']
                if hessian.spacial_dof != recipe.dof:
                    return exit_failure('error: hessian shape is incorrect (wrong DOF!)')
                if not args.do_not_steal:
                    if args.verbose > 0:
                        print('!! {} hessian to perform projection (and geometry)'.format(
                            'replacing' if 'GG' in cf.derivatives else 'adding'))
                    cf.derivatives['GG'] = hessian
                    recipe.geometry = args.hessian.molecule
            except (PropertyNotDefined, PropertyNotPresent):
                return exit_failure('error: file does not contain any hessian (or it cannot find it)')
        else:
            if 'GG' in cf.derivatives:
                hessian = cf.derivatives['GG']

        if hessian is None:
            return exit_failure('error: projection requested, but no hessian available')

        mwh = derivatives_g.MassWeightedHessian(recipe.geometry, cartesian_hessian=hessian)

        try:
            baking.project_geometrical_derivatives(recipe, cf, mwh, verbosity_level=args.verbose)
        except ValueError as e:
            return exit_failure('error while projecting: {}'.format(str(e)))

    with open(args.output, 'w') as f:
        cf.write(f)


if __name__ == '__main__':
    main()
