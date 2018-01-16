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
{rm}
.MXLRV
{max_it}
.THRENR
{cc_conv}
.THRLEQ
{resp_thr}"""

WAVE_FUNCTION = """
**WAVE FUNCTIONS
*SCF INPUT
.THRESH
{conv}
"""


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
                if not compute_FDx and (basis in ['XDD', 'dDF']):
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
        files_created = 0
        inputs_matching = ''

        dal_files = {}
        if 'frequencies' in self.recipe:
            frequencies = ' '.join(
                '{:.8f}'.format(derivatives_e.convert_frequency_from_string(a)) for a in self.recipe['frequencies'])

            freq_card = dalton.InputCard(parameters=[len(self.recipe['frequencies']), frequencies])
        else:
            frequencies = ''
            freq_card = None

        thclr_card = dalton.InputCard(parameters=[
            '{:.1e}'.format(float(self.recipe['flavor_extra']['response_threshold'])).replace('e', 'D')])

        for fields, level in self.fields_needed_by_recipe:
            counter += 1

            bases = self.recipe.bases(level_min=level)

            # separate in group
            groups = {}
            for b in (b[0] for b in bases):
                if b.order() not in groups:
                    groups[b.order()] = []
                groups[b.order()].append(b)

            # merge energy, eventual order 1 and 2
            m_group = []
            x_merge = [0, 1]
            x_split = []

            if int(self.recipe['flavor_extra']['split_level_3']) != 0:
                x_split.append(3)

            if int(self.recipe['flavor_extra']['split_level_4']) != 0:
                x_split.append(4)

            if self.recipe['method'] == 'CC':
                x_merge = [0, 1, 2]

                if int(self.recipe['flavor_extra']['merge_level_3']) != 0:
                    x_merge.append(3)

                if int(self.recipe['flavor_extra']['merge_level_4']) != 0:
                    x_merge.append(4)

            for i in x_merge:
                if i in groups:
                    m_group.extend(groups[i])
                    del groups[i]

            if m_group:
                groups[0] = m_group

            # split groups
            for i in x_split:
                if i in groups:
                    for j, d in enumerate(groups[i]):
                        groups[i * 100 + j] = [d]
                    del groups[i]

            if 2 in groups and self.recipe['method'] != 'CC':
                if 'GG' in groups[2] and len(groups[2]) > 1:
                    i = groups[2].index('GG')
                    groups[201] = [groups[2][i]]
                    del groups[2][i]

                    if 'G' in groups[0]:  # since hessian compute also gradient, no need to do it here
                        i = groups[0].index('G')
                        del groups[0][i]

            bases_reprs = []
            for _, bases in groups.items():
                bases_reprs.append(tuple([b.representation() for b in bases]))

            bases_reprs.sort(key=lambda x: (len(x[0]), x[0].count('D')))  # sort by complexity

            # dal file, if needed
            for bases_repr in bases_reprs:
                if bases_repr not in dal_files:
                    dal = dalton.Input()

                    if self.recipe['method'] == 'CC':
                        if 'dDFF' in bases_repr:
                            raise BadPreparation('dDFF is not available if CC')

                        # beginning
                        dal.update('**DALTON INPUT\n.RUN WAVE FUNCTIONS')

                        # coupled cluster
                        if self.recipe['method'] != 'CC3':
                            dal['DALTON']['.DIRECT'] = dalton.InputCard()

                        dal.update(CC_WAVE_FUNCTION.format_map({
                            'method': self.recipe['flavor_extra']['CC'],
                            'max_it': self.recipe['flavor_extra']['response_max_it'],
                            'rm': self.recipe['flavor_extra']['response_dim_reduced_space'],
                            'conv': '{:.1e}'.format(float(self.recipe['flavor_extra']['threshold'])).replace('e', 'D'),
                            'cc_conv': '{:.1e}'.format(
                                float(self.recipe['flavor_extra']['cc_threshold'])).replace('e', 'D'),
                            'resp_thr': '{:.1e}'.format(
                                float(self.recipe['flavor_extra']['response_threshold'])).replace('e', 'D')
                        }))

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
                        if any([x in bases_repr for x in ['FF', 'dD']]):
                            dal.update('**WAVE FUNC\n*CCLR\n.DIPOLE\n.STATIC')

                            if 'dD' in bases_repr:
                                dal['WAVE FUNC']['CCLR']['.FREQUE'] = copy.copy(freq_card)

                        # first hyperpolarizability
                        if any([x in bases_repr for x in ['FFF', 'dDF', 'XDD']]):
                            dal.update('**WAVE FUNC\n*CCQR\n.DIPOLE')

                            if 'FFF' in bases_repr:
                                dal['WAVE FUNC']['CCQR']['.STATIC'] = dalton.InputCard()

                            if 'dDF' in bases_repr:
                                dal['WAVE FUNC']['CCQR']['.EOPEFR'] = copy.copy(freq_card)

                            if 'XDD' in bases_repr:
                                dal['WAVE FUNC']['CCQR']['.SHGFRE'] = copy.copy(freq_card)

                        # second hyperpolarizability
                        if any([x in bases_repr for x in ['FFFF', 'dFFD', 'XDDF', 'dDDd', 'XDDD']]):
                            dal.update('**WAVE FUNC\n*CCCR\n.DIPOLE')

                            if 'FFFF' in bases_repr:
                                dal['WAVE FUNC']['CCCR']['.STATIC'] = dalton.InputCard()

                            if 'dFFD' in bases_repr:
                                dal['WAVE FUNC']['CCCR']['.DCKERR'] = copy.copy(freq_card)

                            if 'XDDF' in bases_repr:
                                dal['WAVE FUNC']['CCCR']['.ESHGFR'] = copy.copy(freq_card)

                            if 'dDDd' in bases_repr:
                                dal['WAVE FUNC']['CCCR']['.DFWMFR'] = copy.copy(freq_card)

                            if 'XDDD' in bases_repr:
                                dal['WAVE FUNC']['CCCR']['.THGFRE'] = copy.copy(freq_card)

                    else:
                        if 'dFFD' in bases_repr:
                            raise BadPreparation('dFFD is not available if not CC')

                        dal.update('**DALTON INPUT\n.DIRECT\n.RUN WAVE FUNCTIONS')

                        dal.update(WAVE_FUNCTION.format_map({
                            'conv': '{:.1e}'.format(float(self.recipe['flavor_extra']['threshold'])).replace('e', 'D')
                        }))

                        if self.recipe['method'] == 'HF':
                            dal.update('**WAVE FUNCTIONS\n.HF')
                        else:
                            dal.update('**WAVE FUNCTIONS\n.DFT\n{}'.format(self.recipe['flavor_extra']['XC']))
                            # ???

                        # dipole moment
                        if 'F' in bases_repr:
                            dal.update('**DALTON INPUT\n.RUN PROPERTIES')
                            dal.update('**INTEGRAL\n.DIPLEN\n')

                        # gradient and hessian
                        if 'G' in bases_repr and 'GG' not in bases_repr:
                            dal.update('**DALTON INPUT\n.RUN PROPERTIES')
                            dal.update('**PROPERTIES\n.MOLGRA')

                        if 'GG' in bases_repr:
                            dal.update('**DALTON INPUT\n.RUN PROPERTIES')
                            dal.update('**PROPERTIES\n.VIBANA\n*VIBANA\n.HESPUN\n*RESPON')
                            dal['PROPERTIES']['RESPON']['.THRESH'] = copy.copy(thclr_card)

                        if any([x in bases_repr for x in [
                                'FF', 'dD', 'FFF', 'dDF', 'XDD', 'FFFF', 'dFFD', 'XDDF', 'XDDD', 'dDDd']]):
                            dal.update('**DALTON INPUT\n.RUN RESPONSE')
                            dal.update('**RESPONSE')
                            dal['RESPONSE']['.MAXRM'] = dalton.InputCard(parameters=[
                                self.recipe['flavor_extra']['response_dim_reduced_space']])

                        # polarizability
                        if any([x in bases_repr for x in ['FF', 'dD']]):
                            dal.update('**RESPONSE\n*LINEAR\n.DIPLEN')

                            dal['RESPONSE']['LINEAR']['.THCLR'] = copy.copy(thclr_card)
                            dal['RESPONSE']['LINEAR']['.MAX IT'] = dalton.InputCard(
                                parameters=[self.recipe['flavor_extra']['response_max_it']])
                            dal['RESPONSE']['LINEAR']['.MAXITO'] = dalton.InputCard(
                                parameters=[self.recipe['flavor_extra']['response_max_ito']])

                            fx = frequencies
                            num_frequencies = 0
                            if 'frequencies' in self.recipe:
                                num_frequencies = len(self.recipe['frequencies'])
                            if 'FF' in bases_repr:
                                fx = '0.0 ' + fx
                                num_frequencies += 1

                            dal['RESPONSE']['LINEAR']['.FREQUE'] = dalton.InputCard(parameters=[num_frequencies, fx])

                        # first hyperpolarizability
                        if any([x in bases_repr for x in ['FFF', 'dDF', 'XDD']]):
                            dal.update('**RESPONSE\n*QUADRATIC\n.DIPLEN')

                            dal['RESPONSE']['QUADRATIC']['.THCLR'] = copy.copy(thclr_card)
                            dal['RESPONSE']['QUADRATIC']['.MAX IT'] = dalton.InputCard(
                                parameters=[self.recipe['flavor_extra']['response_max_it']])
                            dal['RESPONSE']['QUADRATIC']['.MAXITO'] = dalton.InputCard(
                                parameters=[self.recipe['flavor_extra']['response_max_ito']])

                            fx = frequencies
                            num_frequencies = 0
                            if 'frequencies' in self.recipe:
                                num_frequencies = len(self.recipe['frequencies'])
                            if 'FFF' in bases_repr:
                                fx = '0.0 ' + fx
                                num_frequencies += 1

                            if any([x in bases_repr for x in ['dDF', 'XDD']]):
                                dal['RESPONSE']['QUADRATIC']['.FREQUE'] = dalton.InputCard(
                                    parameters=[num_frequencies, fx])

                                if 'dDF' in bases_repr:
                                    dal['RESPONSE']['QUADRATIC']['.POCKEL'] = dalton.InputCard()
                                if 'XDD' in bases_repr:
                                    dal['RESPONSE']['QUADRATIC']['.SHG'] = dalton.InputCard()

                        # second hyperpolarizability
                        if any([x in bases_repr for x in ['FFFF', 'dDFF', 'XDDF', 'dDDd', 'XDDD']]):
                            dal.update('**RESPONSE\n*CUBIC\n.DIPLEN\n.GAMALL')

                            dal['RESPONSE']['CUBIC']['.THCLR'] = copy.copy(thclr_card)
                            dal['RESPONSE']['CUBIC']['.MAX IT'] = dalton.InputCard(
                                parameters=[self.recipe['flavor_extra']['response_max_it']])
                            dal['RESPONSE']['CUBIC']['.MAXITO'] = dalton.InputCard(
                                parameters=[self.recipe['flavor_extra']['response_max_ito']])

                            fx = frequencies
                            num_frequencies = 0
                            if 'frequencies' in self.recipe:
                                num_frequencies = len(self.recipe['frequencies'])
                            if 'FFFF' in bases_repr:
                                fx = '0.0 ' + fx
                                num_frequencies += 1

                            if any([x in bases_repr for x in ['dDFF', 'XDDF', 'dDDd', 'XDDD']]):
                                dal['RESPONSE']['CUBIC']['.FREQUE'] = dalton.InputCard(
                                    parameters=[num_frequencies, fx])

                                if 'dDFF' in bases_repr:
                                    dal.update('**RESPONSE\n*CUBIC\n.DC-KERR')

                                if 'XDDF' in bases_repr:
                                    dal.update('**RESPONSE\n*CUBIC\n.DC-SHG')

                                if 'XDDD' in bases_repr:
                                    dal.update('**RESPONSE\n*CUBIC\n.THG')

                                if 'dDDd' in bases_repr:
                                    dal.update('**RESPONSE\n*CUBIC\n.IDRI')

                    dal_path = '{}_{}.dal'.format(self.recipe['flavor_extra']['dal_name'], '_'.join(
                        a if a != '' else 'energy' for a in bases_repr))
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

                for bases_repr in bases_reprs:
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
