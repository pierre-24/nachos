from qcip_tools import chemistry_files, molecule, atom
from qcip_tools.chemistry_files import helpers, PropertyNotDefined


class ProgramCalled:
    """Remind when a program is called

    :param prog_name: program
    :type prog_name: str
    :type line_start: int
    :type line_end: int
    """

    def __init__(self, prog_name, line_start, line_end):

        self.prog_name = prog_name
        self.line_start = line_start
        self.line_end = line_end

    def __repr__(self):
        return 'Program {}: {}:{}'.format(self.prog_name, self.line_start, self.line_end)


class QChemLogFile(
        chemistry_files.ChemistryLogFile, chemistry_files.WithIdentificationMixin, chemistry_files.WithMoleculeMixin):

    #: The identifier
    file_type = 'QCHEM_LOG'
    chunk_title_variable = 'prog_name'

    def __init__(self):
        self.molecule = molecule.Molecule()

    @classmethod
    def attempt_identification(cls, f):
        """A QChem log ... Contains a few "qchem" in the beginning (limit to the 100 first lines)"
        """

        count = 0
        num_of_qchem = 0

        for line in f.readlines():
            if count > 100:
                break
            if 'Q-Chem' in line or 'qchem' in line:
                num_of_qchem += 1
            count += 1

        return num_of_qchem > 5

    def read(self, f):
        """

        :param f: File
        :type f: file
        """

        super().read(f)
        self.chunks.append(ProgramCalled('intro', 1, -1))

        line_geom = -1

        for index, line in enumerate(self.lines):
            if 'General SCF calculation program' in line:
                self.chunks[-1].line_end = index - 1
                self.chunks.append(ProgramCalled('scf', index, -1))
            elif 'Standard Nuclear Orientation (Angstroms)' in line:
                line_geom = index + 2
            elif 'CCMAN2:' in line:
                self.chunks[-1].line_end = index - 1
                self.chunks.append(ProgramCalled('ccman2', index, -1))
            elif 'Orbital Energies (a.u.) and Symmetries' in line:
                self.chunks[-1].line_end = index - 1
                self.chunks.append(ProgramCalled('analysis', index, -1))

        if line_geom == -1:
            raise Exception('geometry not found')

        for line in self.lines[line_geom:]:
            if '---------------' in line:
                break

            inf = line.split()
            self.molecule.insert(atom.Atom(
                symbol=inf[1],
                position=[float(a) for a in inf[2:5]]
            ))


# add to available helpers:
helpers.EXTRA_CHEMISTRY_FILES.append(QChemLogFile)


# add properties
@QChemLogFile.define_property('n:input_electric_field')
def qchem__log__property__input_electric_field(obj, *args, **kwargs):
    """Get the input electric field

    :param obj: object
    :type obj: QChemLogFile
    :rtype: dict
    """

    field = [.0, .0, .0, .0]

    found = obj.search('Cartesian multipole field', into='intro')
    if found == -1:
        return field

    for i in range(3):
        field[i + 1] = float(obj.lines[found + 3 + i][16:].strip())

    return field


@QChemLogFile.define_property('computed_energies')
def qchem__log__property__computed_energies(obj, *args, **kwargs):
    """Get the energies. Returns a dictionary of the energies at different level of approximation.

    :param obj: object
    :type obj: QChemLogFile
    :rtype: dict
    """

    found = obj.search('SCF energy', into='ccman2')
    if found == -1:
        raise PropertyNotDefined('computed_energy')

    energies = {'total': .0}

    for line in obj.lines[found:]:
        if len(line) < 5:
            break

        inf = line.split()

        if inf[1] == 'correlation':
            continue

        if inf[0] == 'SCF':
            inf[0] = 'HF'

        e = float(inf[-1])
        energies[inf[0] if inf[0] != 'SCF' else 'HF'] = e
        energies['total'] = e

    return energies
