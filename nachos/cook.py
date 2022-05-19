"""
Out of the results of calculation, create a h5 file to store them
"""

import os
import argparse

import nachos
from nachos.core import files, cooking, preparing
from nachos import exit_failure

__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'


def is_dir(dirname):
    """Checks if a path is an actual directory"""
    if not os.path.isdir(dirname):
        msg = '{0} is not a directory'.format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname


# program options
def get_arguments_parser():
    arguments_parser = argparse.ArgumentParser(description=__doc__)
    arguments_parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    arguments_parser.add_argument('-V', '--verbose', type=int, help='Level of details (0 to 1)', default=0)
    arguments_parser.add_argument(
        '-r', '--recipe', type=argparse.FileType('r'), help='Recipe file', default='./nachos_recipe.yml')
    arguments_parser.add_argument(
        '-o', '--output', type=str, help='Output h5 file', default='nachos_data.h5')

    arguments_parser.add_argument('directories', nargs='*', type=is_dir, help='directory where to look for QM results')

    arguments_parser.add_argument(
        '--gaussian-logs', action='store_true', help='Use Gaussian LOGs instead of FCHKs (... but why in the world?!?)')

    return arguments_parser


# main
def main():
    args = get_arguments_parser().parse_args()

    if args.verbose >= 1:
        print('This is nachos_cook (v{}) from nachos (v{})\n'.format(__version__, nachos.__version__))

    recipe_directory = os.path.dirname(args.recipe.name)
    recipe = files.Recipe(directory=recipe_directory)

    # read recipe and cook
    try:
        recipe.read(args.recipe)
    except files.BadRecipe as e:
        return exit_failure('error while opening recipe: {}'.format(str(e)))

    cooker = cooking.Cooker(recipe)
    directories = [recipe_directory]

    if len(args.directories) != 0:
        directories = args.directories

    try:
        storage = cooker.cook(directories, verbosity_level=args.verbose, use_gaussian_logs=args.gaussian_logs)
    except cooking.BadCooking as e:
        return exit_failure('error while cooking inputs: {}'.format(str(e)))

    missing_fields, missing_derivatives = storage.check()

    if len(missing_fields) != 0 or len(missing_derivatives) != 0:
        num_missing_fields, num_missing_derivs = 0, 0
        errors = ''
        for f in missing_fields:
            errors += '- Missing field: {}\n'.format(','.join(
                preparing.Preparer.nonzero_fields(f, recipe.geometry, recipe['type'])))
            num_missing_fields += 1
        for f, d in missing_derivatives:
            errors += '- Missing derivative {} for {}\n'.format(d, ','.join(
                preparing.Preparer.nonzero_fields(f, recipe.geometry, recipe['type'])))
            num_missing_derivs += 1
        errors += '-' * 32 + '\n'
        errors += 'Total: {} missing field(s) and {} missing derivative(s)'.format(
            num_missing_fields, num_missing_derivs)
        return exit_failure('Errors:\n{}'.format(errors))
    else:
        storage.write(args.output)


if __name__ == '__main__':
    main()
