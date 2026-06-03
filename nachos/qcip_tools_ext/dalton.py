import numpy

from qcip_tools.chemistry_files import dalton
from qcip_tools.chemistry_files import PropertyNotPresent


@dalton.ArchiveOutput.define_property('n:input_electric_field')
def gaussian__FCHK__get_input_electric_field(obj: dalton.ArchiveOutput, *args, **kwargs) -> numpy.ndarray:
    """
    Get the input electric field, look for `@  The molecule is placed in a static field` in `DALTON.CM`.
    """

    try:
        f = obj.get_file('DALTON.CM')
    except FileNotFoundError:
        raise PropertyNotPresent('n:input_electric_field')

    lines = f.readlines()
    field = numpy.zeros(4)

    found_field = -1
    for i, line in enumerate(lines):
        if '@  The molecule is placed in a static field' in line.decode():
            found_field = i

    if found_field > 0:
        p = {'X': 1, 'Y': 2, 'Z': 3}
        for line in lines[found_field + 3:]:
            line = line.decode()
            if '-----' in line:
                break
            else:
                chunks = line.split()
                field[p[chunks[2][0]]] = float(chunks[1])

    return field
