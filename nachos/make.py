"""
Make a recipe
"""

import argparse

from nachos import exit_failure
from nachos.core import making, files

__version__ = '0.1'
__author__ = 'Pierre Beaujean'
__maintainer__ = 'Pierre Beaujean'
__email__ = 'pierre.beaujean@unamur.be'
__status__ = 'Development'


# program options
def get_arguments_parser():
    arguments_parser = argparse.ArgumentParser(description=__doc__)
    arguments_parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    arguments_parser.add_argument(
        '-N', '--fallback-prompt', action='store_true', help='do not use prompt_toolkit')
    arguments_parser.add_argument(
        '-S', '--stop-when-fail', action='store_true', help='fail when argument input is wrong')

    arguments_parser.add_argument('--flavor', help='flavor for the recipe')
    arguments_parser.add_argument('--type', help='type of differentiation')
    arguments_parser.add_argument('--method', help='computational method')
    arguments_parser.add_argument('--basis-set', help='basis set')
    arguments_parser.add_argument('--geometry', help='geometry of the molecule')
    arguments_parser.add_argument('--differentiation', help='differentiation')
    arguments_parser.add_argument('--frequencies', help='frequencies (if dynamic quantities)')
    arguments_parser.add_argument('--name', help='Name of the files')
    arguments_parser.add_argument('--min-field', help='Minimum field (F_0)')
    arguments_parser.add_argument('--ratio', help='ratio (a)')
    arguments_parser.add_argument('--k-max', help='Maximum k (k_max)')
    arguments_parser.add_argument('--flavor-extra', help='Update the values of flavor extra')

    arguments_parser.add_argument('--XC', help='XC functional (if DFT)')
    arguments_parser.add_argument('--gen-basis', help='gaussian basis function file (if gen)')

    arguments_parser.add_argument(
        '-o', '--output', type=argparse.FileType('w'), help='Output recipe file', default='./nachos_recipe.yml')

    return arguments_parser


# main
def main():
    args = get_arguments_parser().parse_args()

    maker = making.Maker(
        use_fallback_prompt=args.fallback_prompt, raise_when_arg_wrong=args.stop_when_fail)
    try:
        recipe = maker.make(args)
    except KeyboardInterrupt:
        return
    except making.BadMaking as e:
        return exit_failure('error while making recipe: {}'.format(str(e)))

    try:
        recipe.check_data()
    except files.BadRecipe as e:
        return exit_failure('error while making recipe: {}'.format(str(e)))

    recipe.write(args.output)

if __name__ == '__main__':
    main()
