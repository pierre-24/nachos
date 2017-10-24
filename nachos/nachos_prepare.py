"""
Create the input files for the quantum chemistry programs out of a given recipe
"""

import os
import argparse
import shutil

from nachos.core import files, exit_failure, cooking

__version__ = '0.2'
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

    arguments_parser.add_argument(
        '-d', '--directory', action='store', help='output directory', default='.', type=is_dir)
    arguments_parser.add_argument(
        '-r', '--recipe', type=argparse.FileType('r'), help='Recipe file', default='nachos_recipe.yml')
    arguments_parser.add_argument(
        '-c', '--copy-files', action='store_true', help='copy geometry and recipe into destination directory')

    return arguments_parser


# main
def main():
    args = get_arguments_parser().parse_args()

    recipe = files.Recipe()

    # move chdir so that it correspond to the directory where the recipe is
    directory = os.path.abspath(args.directory)
    chdir = os.path.dirname(args.recipe.name)
    os.chdir(chdir)

    # read recipe and prepare!
    try:
        recipe.read(args.recipe)
    except files.BadRecipe as e:
        return exit_failure('error while opening recipe: {}'.format(str(e)))

    preparer = cooking.Preparer(recipe, directory)

    try:
        n = preparer.prepare()
    except cooking.BadPreparation as e:
        return exit_failure('error while cooking inputs: {}'.format(str(e)))

    if args.copy_files:
        shutil.copy(recipe['geometry'], directory)
        recipe['geometry'] = os.path.basename(recipe['geometry'])

        if 'gen_basis' in recipe['flavor_extra']:
            shutil.copy(recipe['flavor_extra']['gen_basis'], directory)
            recipe['flavor_extra']['gen_basis'] = os.path.basename(recipe['flavor_extra']['gen_basis'])

        os.chdir(directory)
        with open('nachos_recipe.yml', 'w') as f:
            recipe.write(f)
        os.chdir(chdir)

    print('cooked {} files'.format(n))

if __name__ == '__main__':
    main()
