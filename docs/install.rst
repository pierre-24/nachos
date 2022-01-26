=================
Installing nachos
=================

Normal install
--------------

To install the latest version of nachos:

.. code-block:: bash

  # install qcip_tools (base library):
  pip3 install --user --upgrade git+ssh://git@gitlab.unamur.be/chimie/lct/qcip_tools.git@release-v0.6.1
  # install nachos:
  pip3 install --user --upgrade git+ssh://git@gitlab.unamur.be/chimie/lct/nachos.git

Add ``@release-vXX`` at the end of the last line to fetch a given version (listed `in the README <https://gitlab.unamur.be/chimie/lct/nachos/blob/master/README.md>`_).


Note that ``--user`` allow you to install the package without being superuser (see `here <https://pip.pypa.io/en/stable/user_guide/#user-installs>`_).
You will probably need to add ``$HOME/.local/bin`` to ``$PATH`` for this to work:

.. code-block:: bash

  echo 'PATH=$HOME/.local/bin:$PATH' >> ~/.bashrc

On the other hand, you can install it in a *virtualenv* (see below).


(*optional*) Patching Dalton
----------------------------

By default, it is possible to perform numerical differentiation with Dalton, but the following patch will improve different stuffs:

+ Increase values of some constants, so that more frequencies and responses functions can be computed in the same input ;
+ Outputs responses functions in a better place (``DALTON.PROP`` in the archive), with more digits (important for the accuracy) ;
+ Allow to compute numerical differentiation of gamma (because it is otherwise not possible, since only the components that participate to :math:`\gamma_{||}` are computed).

You can download the patch here: :download:`dalton.patch` :

.. code-block:: bash

  wget --http-user=doc --http-passwd=doc http://perso.unamur.be/~pbeaujea/nachos/html/_downloads/dalton.patch

To apply this patch you need to recompile Dalton (so you need git, a fortran compiler, and eventual BLAS/LAPACK/MKL).
The following  commands allow you to `get the sources <https://gitlab.com/dalton/dalton>`_, and check if everything is ok:

.. code-block:: bash

  git clone --recursive https://gitlab.com/dalton/dalton.git
  cd dalton/
  git checkout -b special_version origin/release/2016
  # use "origin/master" if you want the latest (maybe buggy) version
  git apply --check /path/to/dalton.patch

If it is the case (no ``error`` in the output of the previous command), you can apply it and compile a new version of Dalton:

.. code-block:: bash

  git am --signoff < /path/to/dalton.patch
  #  You may get "trailing whitespace" error, but it's ok
  ./setup --prefix=/path/to/destination # [--int64]
  cd build
  make
  make install

Installation for contributors
-----------------------------

To contribute to the project, you need to clone the repository:

+ Clone it: ``git clone git@git.pierrebeaujean.net:pierre/nachos.git``.
+ Install pip-tools: ``pip3 install pip-tools``
+ Install virtualenv ``python3 -m venv venv; source venv/bin/activate``
+ Install dependencies: ``make init``.
+ Don't forget to create a separate branch to implement your changes (see `the contribution part <contributing.html>`_).

You can launch the tests series with ``make test``