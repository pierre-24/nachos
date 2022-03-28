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


@gaussian.Output.define_property('computed_energies')
def gaussian__Output__get_computed_energies(obj, *args, **kwargs):
    """Get the energies... Actually, only
    + HF and DFT energy in link 502 (normal SCF, not QC which is link 508)
    + MP2 energy in link 906

    :param obj: object
    :type obj: qcip_tools.chemistry_files.gaussian.Output
    :rtype: dict
    """

    energies = {}

    # fetch HF or DFT energy
    n = obj.search('SCF Done:', into=502)
    if n > 0:
        chunks = obj.lines[n][12:].split()
        e = float(chunks[2])
        if 'HF' in chunks[0]:
            energies['HF'] = e
        else:
            energies['SCF/DFT'] = e

        energies['total'] = e

    # fetch MP2 energies
    n = obj.search('EUMP2 =', into=906)
    if n > 0:
        energies['MP2'] = float(obj.lines[n].split()[-1])
        energies['total'] = energies['MP2']

    return energies


@gaussian.Output.define_property('n:spin_components_e2')
def gaussian__Output__get_spin_component_e2(obj, *args, **kwargs):
    """Get the spin components of E(2) (link 906)

    :param obj: object
    :type obj: qcip_tools.chemistry_files.gaussian.Output
    :rtype: dict
    """

    n = obj.search('Spin components of T(2) and E(2)', into=906)
    if n > 0:
        return {
            'aa': float(obj.lines[n + 1][53:-1].replace('D', 'E')),
            'ab': float(obj.lines[n + 2][53:-1].replace('D', 'E')),
            'bb': float(obj.lines[n + 3][53:-1].replace('D', 'E'))  # if closed shell, aa = bb
        }
    else:
        raise PropertyNotPresent('n:spin_components_e2')
