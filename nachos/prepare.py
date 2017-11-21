"""
Create the input files for the quantum chemistry programs out of a given recipe
"""

import argparse
import os
import shutil

from nachos import is_dir, exit_failure
from nachos.core import files, preparing

__version__ = '0.2'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'


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
    directory = os.path.abspath(args.directory)
    recipe_directory = os.path.dirname(args.recipe.name)

    recipe = files.Recipe(directory=recipe_directory)

    # read recipe and prepare!
    try:
        recipe.read(args.recipe)
    except files.BadRecipe as e:
        return exit_failure('error while opening recipe: {}'.format(str(e)))

    preparer = preparing.Preparer(recipe, directory)

    try:
        n = preparer.prepare()
    except preparing.BadPreparation as e:
        return exit_failure('error while cooking inputs: {}'.format(str(e)))

    if args.copy_files:
        shutil.copy(os.path.join(recipe_directory, recipe['geometry']), directory)
        recipe['geometry'] = os.path.basename(os.path.join(recipe_directory, recipe['geometry']))

        if 'gen_basis' in recipe['flavor_extra']:
            shutil.copy(os.path.join(recipe_directory, recipe['flavor_extra']['gen_basis']), directory)
            recipe['flavor_extra']['gen_basis'] = os.path.basename(
                os.path.join(recipe_directory, recipe['flavor_extra']['gen_basis']))

        with open(os.path.join(directory, 'nachos_recipe.yml'), 'w') as f:
            recipe.write(f)

    print('cooked {} files'.format(n))

if __name__ == '__main__':
    main()
