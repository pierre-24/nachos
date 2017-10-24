import sys

CONFIG = {
    'gaussian': {
        'types': ['G', 'F'],
        'methods': [
            ('HF', {'G': 2, 'F': 3}),
            ('DFT', {'G': 2, 'F': 3}),
            ('MP2', {'G': 2, 'F': 2}),  # ... but only static polarizability
            ('MP3', {'G': 1, 'F': 1}),
            ('MP4D', {'G': 1, 'F': 1}),
            ('MP4DQ', {'G': 1, 'F': 1}),
            ('MP4SDQ', {'G': 1, 'F': 1}),
            ('MP4', {'G': 1, 'F': 1}),
            ('CCSD', {'G': 1, 'F': 1}),
            ('CCSD(T)', {'G': 1, 'F': 1}),
        ],
        'bases': ['energy', 'G', 'GG', 'F', 'FF', 'FD', 'FFF', 'FDF', 'FDD'],
        'default_for_extra_fields': {
            'memory': '1Gb',
            'procs': 1,
            'convergence': 11,
            'cc_convergence': 11,
            'cphf_convergence': 10,
            'max_cycles': 600,
            'extra_keywords': '',
            'extra_sections': [],
            'gen_basis': '',
            'vshift': 1000,
            'XC': ''
        }
    },
    'dalton': {
        'types': ['G'],
        'methods': [
            ('CCS', {'G': 1, 'F': 4}),
            ('CC2', {'G': 1, 'F': 4}),
            ('CCSD', {'G': 1, 'F': 4}),
            ('CC3', {'G': 1, 'F': 4}),
        ],
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
