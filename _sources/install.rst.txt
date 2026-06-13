=================
Installing nachos
=================

Normal install
--------------

To install the latest version of nachos:

.. code-block:: bash

  pip3 install --user --upgrade git+https://github.com/pierre-24/nachos.git

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

You can find the patch and instructions `there <https://pierre-24.github.io/qcip_tools/install.html#optional-patching-dalton>`_.