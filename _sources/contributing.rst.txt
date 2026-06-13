============
Contributing
============

Installation for contributors
-----------------------------

To contribute to the project,

+ `Fork it <https://docs.github.com/en/get-started/quickstart/fork-a-repo>`_.
+ Clone your fork: ``git clone git@github.com:<USERNAME>/nachos.git``.
+ Go in it: ``cd nachos``
+ Set up a virtualenv ``python3 -m venv venv; source venv/bin/activate``
+ Install dependencies: ``make install-dev``.
+ Add upstream: ``git remote add upstream https://github.com/pierre-24/nachos.git``
+ Don't forget to create a separate branch to implement your changes: ``git checkout -b my_branch upstream/dev``.

You can launch the tests series with ``make test``

You may want to take a look into `qcip_tools <https://github.com/pierre-24/qcip_tools>`_, the base library behind nachos.

Design rules
------------

+ The code is written in Python 3, and follows the (in)famous `PEP-8 <http://legacy.python.org/dev/peps/pep-0008/>`_. You can check it by running ``make lint``, which launch the ``flake`` utility.
+ Codes and comments are written in english.
+ The code is documented using docstrings and Sphinx. The docstrings must contains the basic description of the function, as well as a description of the parameters.
+ The code is tested. You can launch the test series by using ``make test``.
  Every functionality should be provided with at least one unit test.
  Every script should be provided with at least one unit test.
  You may need test files to do so, but try to make them small (say, don't use d-aug-cc-pVDZ while STO-3G could do the job).
+ The package is documented. You can generate this documentation by using ``make doc``. Non-basic stuffs should be explained in this documentation. Don't forget to cite some articles or website if needed.

Workflow
--------

Adapted from the (in)famous `Git flow <http://nvie.com/posts/a-successful-git-branching-model/>`_.

+ Development is made in ``dev`` branch.
+ Functionalities are added through pull requests (PR) to the ``dev`` branch. Do not work in ``dev`` directly, but create a new branch (``git checkout -b my_branch upstream/dev``).
+ Theses pull requests should be unitary, and include unit test(s) and documentation if needed. The test suite must succeed for the merge request to be accepted.
+ The pull requests will be reviewed before acceptance.
+ At some (random) points, a new version will appear, with a tag of the form ``vXX``.