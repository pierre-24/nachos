import copy
import math
import os

import numpy
from qcip_tools import derivatives, derivatives_e, quantities, numerical_differentiation
from qcip_tools.chemistry_files import gaussian, dalton

from nachos.core import compute_numerical_derivative_of_tensor


def fields_needed_by_recipe(recipe):
    """

    :param recipe: recipe
    :type recipe: nachos.core.files.Recipe
    :rtype: list
    """

    diffs = recipe.maximum_derivatives()

    fields_needed = [
        [0] * (3 if recipe['type'] == 'F' else recipe.dof)  # with zero field
    ]

    fields_needed_with_level = [
        (fields_needed[0], 1)
    ]

    def collect_fields(fields, *args, **kwargs):
        f = list(int(a) for a in fields)
        if f not in fields_needed:
            fields_needed.append(f)
            fields_needed_with_level.append((f, kwargs['level']))
        return .0

    e = derivatives.Derivative('')

    for d in diffs:
        compute_numerical_derivative_of_tensor(
            recipe, e, d, collect_fields, dry_run=True, level=d.order())

    return fields_needed_with_level


class BadPreparation(Exception):
    pass


CC_WAVE_FUNCTION = """
**WAVE FUNCTIONS
.CC
*SCF INPUT
.THRESH
{conv}
*CC INPUT
.{method}
.MAX IT
{max_it}
.MAXRED
{max_it}
.MXLRV
{max_it}
.THRENR
{cc_conv}
.THRLEQ
{conv}"""


class Preparer:
    """Prepare the input files

    :param recipe: a recipe
    :type recipe: nachos.core.files.Recipe
    """

    def __init__(self, recipe, directory='.'):
        self.recipe = recipe

        if not os.path.isdir(directory):
            raise BadPreparation('{} is not a directory'.format(directory))

        self.directory = directory
        self.fields_needed_by_recipe = fields_needed_by_recipe(self.recipe)

    def prepare(self):
        """Create the different input files in the directory"""

        return getattr(self, 'prepare_{}_inputs'.format(self.recipe['flavor']))()

    def prepare_gaussian_inputs(self):
        """Create inputs for gaussian
        """

        base_m = False
        counter = 0
        files_created = 0

        # try to open custom basis set if any
        gen_basis_set = []
        if self.recipe['basis_set'] == 'gen':
            path = os.path.join(self.recipe.directory, self.recipe['flavor_extra']['gen_basis'])
            if not os.path.exists(path):
                raise BadPreparation('gen basis file {} cannot be opened'.format(path))

            gbs = gaussian.BasisSet()
            try:
                with open(path) as f:
                    gbs.read(f)
            except gaussian.BasisSetFormatError as e:
                raise BadPreparation('error while opening custom basis set ({}): {}'.format(path, str(e)))

            try:
                gen_basis_set = gbs.to_string(for_molecule=self.recipe.geometry).splitlines()[1:]
            except Exception as e:
                raise BadPreparation('error while using custom basis set ({}) : {}'.format(path, e))

        # try to deal with extra sections if any
        extra_sections = None
        if 'extra_sections' in self.recipe['flavor_extra'] and self.recipe['flavor_extra']['extra_sections'] != '':
            path = os.path.join(self.recipe.directory, self.recipe['flavor_extra']['extra_sections'])
            if not os.path.exists(path):
                raise BadPreparation('extra section file {} cannot be opened'.format(path))

            with open(path) as f:
                content = f.readlines()

            extra_sections = []
            current_section = []
            for l in content:
                c = l.strip()
                if c == '':
                    extra_sections.append(current_section)
                    current_section = []
                else:
                    current_section.append(c)

            if current_section:
                extra_sections.append(current_section)

        for fields, level in self.fields_needed_by_recipe:
            counter += 1

            compute_polar = False
            compute_polar_with_freq = False
            compute_G = False
            compute_GG = False
            compute_FDx = False

            bases = self.recipe.bases(level_min=level)

            for b, l in bases:
                basis = b.representation()
                if not compute_G and basis == 'G':
                    compute_G = True
                if not compute_GG and basis == 'GG':
                    compute_GG = True
                if not compute_polar and 'FF' in basis:
                    compute_polar = True
                if not compute_polar_with_freq and 'D' in basis:
                    compute_polar = True
                    compute_polar_with_freq = True
                if not compute_FDx and (basis in ['FDD', 'FDF']):
                    compute_FDx = True

            compute_polar_and_G = compute_polar and (compute_G or compute_GG)

            fi = gaussian.Input()
            real_fields = numerical_differentiation.real_fields(fields, self.recipe['min_field'], self.recipe['ratio'])

            if self.recipe['type'] == 'G':
                fi.molecule = Preparer.deform_geometry(self.recipe.geometry, real_fields)
            else:
                fi.molecule = self.recipe.geometry

            if not base_m and fields == [0] * len(fields):
                fi.title = 'base'
                base_m = True
            else:
                fi.title = 'field({})='.format(level) + \
                    ','.join(Preparer.nonzero_fields(fields, self.recipe.geometry, self.recipe['type']))

            fi.options['nprocshared'] = self.recipe['flavor_extra']['procs']
            fi.options['mem'] = self.recipe['flavor_extra']['memory']
            fi.options['chk'] = 'xxx'

            # input card
            input_card = [
                '#P {} {} nosym{}'.format(
                    self.recipe['method'] if self.recipe['method'] != 'DFT' else self.recipe['flavor_extra']['XC'],
                    self.recipe['basis_set'],
                    ' field=read' if self.recipe['type'] == 'F' else ''),
                'scf=(Conver={c},NoVarAcc,MaxCyc={m},vshift={v}) IOP(9/6={m},9/9={cc})'.format_map(
                    {
                        'c': self.recipe['flavor_extra']['convergence'],
                        'cc': self.recipe['flavor_extra']['cc_convergence'],
                        'm': self.recipe['flavor_extra']['max_cycles'],
                        'v': self.recipe['flavor_extra']['vshift']
                    })
            ]

            if self.recipe['flavor_extra']['extra_keywords']:
                input_card.extend([self.recipe['flavor_extra']['extra_keywords']])

            fi.input_card = input_card

            # other blocks
            if self.recipe['basis_set'] == 'gen':
                fi.other_blocks.append(gen_basis_set)

            if self.recipe['type'] == 'F':
                fi.other_blocks.append(['\t'.join(['{: .10f}'.format(a) for a in real_fields])])

            if extra_sections:
                fi.other_blocks.extend(extra_sections)

            # write files
            if compute_polar:
                extra_line = 'polar{} cphf=(conver={}{})'.format(
                    '=dcshg' if compute_FDx else '',
                    self.recipe['flavor_extra']['cphf_convergence'],
                    ',rdfreq' if compute_polar_with_freq else ''
                )
                fi.input_card.append(extra_line)

                if compute_polar_with_freq:
                    fi.other_blocks.insert(0, [
                        '{}'.format(derivatives_e.convert_frequency_from_string(a)) for a in
                        self.recipe['frequencies']])
                with open('{}/{}_{:04d}{}.com'.format(
                        self.directory, self.recipe['name'], counter, 'a' if compute_polar_and_G else ''), 'w') as f:
                    fi.write(f)
                    files_created += 1
                fi.input_card.pop(-1)
                fi.other_blocks.pop(0)

            if compute_polar_and_G or not compute_polar:
                if compute_GG:
                    fi.input_card.append('freq cphf=(conver={})'.format(
                        self.recipe['flavor_extra']['cphf_convergence']))
                elif compute_G:
                    fi.input_card.append('force')

                with open('{}/{}_{:04d}{}.com'.format(
                        self.directory, self.recipe['name'], counter, 'b' if compute_polar_and_G else ''), 'w') as f:
                    fi.write(f)
                    files_created += 1

        return files_created

    def prepare_dalton_inputs(self):
        """Create inputs for dalton. Note that it assume geometrical derivatives for the moment
        """

        if self.recipe['type'] != 'G':
            raise BadPreparation('Dalton only works for G!')

        base_m = False
        counter = 0
        dal_counter = 0
        files_created = 0
        inputs_matching = ''

        dal_files = {}
        frequencies = ' '.join(
            '{:.8f}'.format(derivatives_e.convert_frequency_from_string(a)) for a in self.recipe['frequencies'])

        freq_card = dalton.InputCard(parameters=[len(self.recipe['frequencies']), frequencies])

        for fields, level in self.fields_needed_by_recipe:
            counter += 1

            bases = self.recipe.bases(level_min=level)
            bases_repr = tuple([b[0].representation() for b in bases])

            # dal file, if needed
            if bases_repr not in dal_files:
                dal_counter += 1
                dal = dalton.Input()

                # begining
                dal.update('**DALTON INPUT\n.RUN WAVE FUNCTIONS')

                # coupled cluster
                if self.recipe['method'] != 'CC3':
                    dal['DALTON']['.DIRECT'] = dalton.InputCard()

                dal.update(CC_WAVE_FUNCTION.format_map({
                    'method': self.recipe['method'],
                    'max_it': self.recipe['flavor_extra']['max_iteration'],
                    'conv': '{:.1e}'.format(self.recipe['flavor_extra']['threshold']).replace('e', 'D'),
                    'cc_conv': '{:.1e}'.format(self.recipe['flavor_extra']['cc_threshold']).replace('e', 'D')}))

                # dipole moment
                if 'F' in bases_repr:
                    dal.update('**INTEGRAL\n.DIPLEN\n')
                    dal.update('**WAVE FUNCTION\n*CCFOP\n.DIPMOM')

                    if self.recipe['method'] == 'CC3':
                        dal['WAVE FUNC']['CCFOP']['.NONREL'] = dalton.InputCard()

                # gradient
                if 'G' in bases_repr:
                    dal.update('**INTEGRAL\n.DEROVL\n.DERHAM\n')
                    dal.update('**WAVE FUNCTION\n*DERIVATIVES')

                # polarizability
                if any([x in bases_repr for x in ['FF', 'FD']]):
                    dal.update('**WAVE FUNC\n*CCLR\n.DIPOLE\n.STATIC')

                    if 'FD' in bases_repr:
                        dal['WAVE FUNC']['CCLR']['.FREQUE'] = copy.copy(freq_card)

                # first hyperpolarizability
                if any([x in bases_repr for x in ['FFF', 'FDF', 'FDD']]):
                    dal.update('**WAVE FUNC\n*CCQR\n.DIPOLE\n.STATIC')

                    if 'FDF' in bases_repr:
                        dal['WAVE FUNC']['CCQR']['.EOPEFR'] = copy.copy(freq_card)

                    if 'FDD' in bases_repr:
                        dal['WAVE FUNC']['CCQR']['.SHGFRE'] = copy.copy(freq_card)

                # second hyperpolarizability
                if any([x in bases_repr for x in ['FFFF', 'FDFF', 'FDDF', 'FDDd', 'FDDD']]):
                    dal.update('**WAVE FUNC\n*CCCR\n.DIPOLE\n.STATIC')

                    if 'FDFF' in bases_repr:
                        dal['WAVE FUNC']['CCCR']['.DCKERR'] = copy.copy(freq_card)

                    if 'FDDF' in bases_repr:
                        dal['WAVE FUNC']['CCCR']['.ESHGFR'] = copy.copy(freq_card)

                    if 'FDDd' in bases_repr:
                        dal['WAVE FUNC']['CCCR']['.DFWMFR'] = copy.copy(freq_card)

                    if 'FDDD' in bases_repr:
                        dal['WAVE FUNC']['CCCR']['.THGFRE'] = copy.copy(freq_card)

                dal_path = '{}_{:04d}.dal'.format(self.recipe['flavor_extra']['dal_name'], dal_counter)
                with open('{}/{}'.format(self.directory, dal_path), 'w') as f:
                    dal.write(f)

                dal_files[bases_repr] = dal_path

            # molecule file
            fi = dalton.MoleculeInput()
            real_fields = numerical_differentiation.real_fields(fields, self.recipe['min_field'], self.recipe['ratio'])
            fi.molecule = Preparer.deform_geometry(self.recipe.geometry, real_fields)

            if not base_m and fields == [0] * len(fields):
                fi.title = 'base'
                base_m = True
            else:
                fi.title = 'field({})='.format(level) + \
                           ','.join(Preparer.nonzero_fields(fields, self.recipe.geometry, self.recipe['type']))

            fi.basis_set = self.recipe['basis_set']

            mol_path = '{}_{:04d}.mol'.format(self.recipe['name'], counter)
            with open('{}/{}'.format(self.directory, mol_path), 'w') as f:
                fi.write(f, nosym=True)
                inputs_matching += '{} {}\n'.format(dal_files[bases_repr], mol_path)
                files_created += 1

        with open('{}/inputs_matching.txt'.format(self.directory), 'w') as f:
            f.write(inputs_matching)

        return files_created

    @staticmethod
    def deform_geometry(geometry, real_fields, geometry_in_angstrom=True):
        """Create an input for gaussian

        :param real_fields: Real differentiation field
        :type real_fields: list
        :param geometry: geometry do deform
        :type geometry: qcip_tools.molecule.Molecule
        :param geometry_in_angstrom: indicate wheter the geometry is given in Angstrom or not
            (because the field is obviously given in atomic units)
        :type geometry_in_angstrom: bool
        :rtype: qcip_tools.molecule.Molecule
        """

        deformed = copy.deepcopy(geometry)

        for index, atom in enumerate(deformed):
            atom.position += numpy.array(
                [real_fields[index * 3 + i] * (
                    quantities.AuToAngstrom if geometry_in_angstrom else 1.) for i in range(3)]
            )

        return deformed

    @staticmethod
    def nonzero_fields(fields, geometry, t):
        return [
            '{}({:+g}{})'.format(
                '{}{}'.format(
                    geometry[int(math.floor(i / 3))].symbol, int(math.floor(i / 3) + 1))
                if t == 'G' else 'F',
                e,
                derivatives.COORDINATES[i % 3]
            ) for i, e in enumerate(fields) if e != 0]
