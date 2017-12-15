=================
Installing nachos
=================

Normal install
--------------

To install nachos for a day-to-day use:

.. code-block:: bash

  # install qcip_tools (base library):
  pip3 install --user --upgrade git+ssh://git@git.pierrebeaujean.net/pierre/qcip_tools.git@dev
  # install the scripts:
  pip3 install --user --upgrade git+ssh://git@git.pierrebeaujean.net/pierre/nachos.git

Note that ``--user`` allow you to install the package without being superuser (see `here <https://pip.pypa.io/en/stable/user_guide/#user-installs>`_).
You will probably need to add ``$HOME/.local/bin`` to ``$PATH`` for this to work:

.. code-block:: bash

  echo 'PATH=$HOME/.local/bin:$PATH' >> ~/.bashrc

On the other hand, you can install it in a *virtualenv* (see below).

Installation for contributors
-----------------------------

To contribute to the project, you need to clone the repository:

+ Clone it: ``git clone git@git.pierrebeaujean.net:pierre/nachos.git``.
+ Create virtualenv and activate it:

.. code-block:: bash

  virtualenv venv --python=python3
  # activate virtualenv (you need to do that every time)
  source venv/bin/activate

+ Install (dev) dependencies : ``make install-dependencies-dev``.
+ Finally, "install" the pakage: ``pip install --editable .``
+ Don't forget to create a separate branch to implement your changes (see `the contribution part <contributing.html>`_)

You can launch the tests series with ``python3 setup.py test``