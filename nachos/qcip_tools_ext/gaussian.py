import numpy

from qcip_tools.chemistry_files import gaussian
from qcip_tools.chemistry_files import PropertyNotPresent


@gaussian.FCHK.define_property('n:input_electric_field')
def gaussian__FCHK__get_input_electric_field(obj, *args, **kwargs):
    """Get the input electric field

    :param obj: object
    :type obj: qcip_tools.chemistry_files.gaussian.FCHK
    :rtype: dict
    """

    if 'External E-field' in obj:
        field = obj.get('External E-field')
        return field
    else:
        raise PropertyNotPresent('n:input_electric_field')


@gaussian.Output.define_property('n:input_electric_field')
def gaussian__Output__get_input_electric_field(obj, *args, **kwargs):
    """Get the input electric field from link 301

    :param obj: object
    :type obj: qcip_tools.chemistry_files.gaussian.Output
    :rtype: dict
    """

    field = numpy.zeros(4)

    n = obj.search('The following finite field(s) will be applied:', into=301)
    if n > 0:
        field[1:] = list(float(x) for x in obj.lines[n + 1][25:].split())

    return field


def _find_in_links(obj, q, links):
    found = filter(lambda n: n > 0, (obj.search(q, into=link) for link in links))
    try:
        return next(found)
    except StopIteration:
        return -1


@gaussian.Output.define_property('computed_energies')
def gaussian__Output__get_computed_energies(obj, *args, **kwargs):
    """Get the energies... Actually, "only"
    + HF and DFT energy (in link 502/503/506/508)
    + MP2 energy (in link 804/903/905/906)
    + MP3, CCSD and CCSD(T) energy (in link 913)

    :param obj: object
    :type obj: qcip_tools.chemistry_files.gaussian.Output
    :rtype: dict
    """
    energies = {}
    modified_gaussian = False  # pristine Gaussian 16

    # fetch HF or DFT energy
    n = _find_in_links(obj, 'SCF Done:', [502, 503, 506, 508])
    if n > 0:
        chunks = obj.lines[n][12:].split()
        e = float(chunks[2])
        if (len(chunks[2]) - chunks[2].find('.')) != 9:
            modified_gaussian = True
        else:  # try to fetch more decimals above: "E = xxx" contains eleven of them
            for line in reversed(obj.lines[:n]):
                if 'Delta-E=' in line:
                    e = float(line.split()[1])
                    break
                elif 'Cycle' in line:
                    break

        if 'HF' in chunks[0]:
            energies['HF'] = e
        else:
            energies['SCF/DFT'] = e

        energies['total'] = e

    # fetch MP2 energies
    n = _find_in_links(obj, 'EUMP2 =', [804, 903, 905, 906])
    if n > 0:
        chunk = obj.lines[n].split()[-1]
        if not modified_gaussian:
            chunk = chunk.replace('D', 'E')

        energies['MP2'] = float(chunk)
        energies['total'] = energies['MP2']

    # fetch MP3 energies
    n = _find_in_links(obj, 'EUMP3=', [913])
    if n > 0:
        chunk = obj.lines[n].split()[-1]
        if not modified_gaussian:
            chunk = chunk.replace('D', 'E')
        energies['MP3'] = float(chunk)

    # fetch CCSD energies
    n = obj.search('Wavefunction amplitudes converged.', into=913)
    if n > 0:
        chunk = obj.lines[n - 3][:-16].split()[-1]
        if not modified_gaussian:
            chunk = chunk.replace('D', 'E')
        energies['CCSD'] = float(chunk)
        energies['total'] = energies['CCSD']

    # fetch CCSD(T) energies
    n = obj.search('CCSD(T)=', into=913)
    if n > 0:
        chunk = obj.lines[n].split()[-1]
        if not modified_gaussian:
            chunk = chunk.replace('D', 'E')
        energies['CCSD(T)'] = float(chunk)
        energies['total'] = energies['CCSD(T)']

    return energies


@gaussian.Output.define_property('n:spin_components_e2')
def gaussian__Output__get_spin_component_e2(obj, *args, **kwargs):
    """Get the spin components of E(2) (link 906)

    :param obj: object
    :type obj: qcip_tools.chemistry_files.gaussian.Output
    :rtype: dict
    """

    n = _find_in_links(obj, 'Spin components of T(2) and E(2)', [804, 903, 905, 906])
    if n > 0:
        return {
            'aa': float(obj.lines[n + 1][53:-1].replace('D', 'E')),
            'ab': float(obj.lines[n + 2][53:-1].replace('D', 'E')),
            'bb': float(obj.lines[n + 3][53:-1].replace('D', 'E'))  # if closed shell, aa = bb
        }
    else:
        raise PropertyNotPresent('n:spin_components_e2')