import sys

CONFIG = {
    'gaussian': {
        'types': ['geometric', 'electric'],
        'methods': ['HF', 'MP2', 'MP3', 'MP4D', 'MP4DQ', 'MP4SDQ', 'MP4', 'CCSD', 'CCSD(T)', 'DFT'],
        'bases': ['energy', 'G', 'GG', 'F', 'FF', 'FD', 'FFF', 'FDF', 'FDD'],
        'default_for_extra_fields': {
            'memory': '1Gb',
            'procs': 1,
            'convergence': 11,
            'cc_convergence': 0,
            'max_cycles': 600,
            'extra_keywords': '',
            'extra_sections': [],
            'gen_basis': '',
            'vshift': 1000,
        }
    },
    'dalton': {
        'types': ['geometric'],
        'methods': ['CCS', 'CC2', 'CCSD', 'CC3'],
        'bases': ['energy', 'G', 'F', 'FF', 'FD', 'FFF', 'FDF', 'FDD', 'FFFF', 'FDFF', 'FDDF', 'FDDd', 'FDDD'],
        'default_for_extra_fields': {
            'max_iteration': 2500,
            'threshold': 1e-6,
            'cc_threshold': 1e-11,
        }
    }
}


def exit_failure(msg, status=1):
    """Write a message in stderr and exits

    :param msg: the msg
    :type msg: str
    :param status: exit status (!=0)
    :type status: int
    """

    sys.stderr.write(msg)
    sys.stderr.write('\n')
    return sys.exit(status)
