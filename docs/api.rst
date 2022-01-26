==================
Code documentation
==================

Documentation for the different part of the "core" nachos package

Core (``nachos.core``)
----------------------

.. automodule:: nachos.core
    :members:


Files (``nachos.core.files``)
-----------------------------

.. automodule:: nachos.core.files
    :members:


Making (``nachos.core.making``)
-------------------------------

.. automodule:: nachos.core.making
    :members:


Preparing (``nachos.core.preparing``)
-------------------------------------

.. automodule:: nachos.core.preparing
    :members:

Cooking (``nachos.core.cooking``)
---------------------------------

.. automodule:: nachos.core.cooking
    :members:


Baking (``nachos.core.baking``)
-------------------------------

.. automodule:: nachos.core.baking
    :members:


Shaking (``nachos.core.shaking``)
---------------------------------

.. note::

    .. math::

        \newcommand{\tdiff}[2]{\left(\frac{\partial #1}{ \partial #2}\right)}
        \newcommand{\lb}[2]{\lambda^{\pm #1}_{#2}}

    The (pv) contributions are written listing the properties, then the anharmonicity orders, so :math:`[\mu^2]^{2,0}` is written ``F_F___2_0``.
    The ZPVA contributions are written the same way, except that it is only ``XDD__0_1`` for :math:`\Delta\beta^{0,1}`.

Here is the list of contributions and corresponding functions:

.. list-table::
   :header-rows: 1

   * - Property
     - pv
     - Function
   * - Polarizability
     - :math:`[\mu^2]^{0,0}`
     - `_compute_F_F__0_0_component <#nachos.core.shaking.Shaker._compute_F_F__0_0_component>`_
   * -
     - :math:`[\mu^2]^{1,1}`
     - `_compute_F_F__1_1_component <#nachos.core.shaking.Shaker._compute_F_F__1_1_component>`_
   * -
     - :math:`[\mu^2]^{2,0}`
     - `_compute_F_F__2_0_component <#nachos.core.shaking.Shaker._compute_F_F__2_0_component>`_
   * -
     - :math:`[\mu^2]^{0,2}`
     - `_compute_F_F__0_2_component <#nachos.core.shaking.Shaker._compute_F_F__0_2_component>`_
   * - First hyperpolarizability
     - :math:`[\mu\alpha]^{0,0}`
     - `_compute_F_FF__0_0_component <#nachos.core.shaking.Shaker._compute_F_FF__0_0_component>`_
   * -
     - :math:`[\mu\alpha]^{1,1}`
     - `_compute_F_FF__1_1_component <#nachos.core.shaking.Shaker._compute_F_FF__1_1_component>`_
   * -
     - :math:`[\mu\alpha]^{2,0}`
     - `_compute_F_FF__2_0_component <#nachos.core.shaking.Shaker._compute_F_FF__2_0_component>`_
   * -
     - :math:`[\mu\alpha]^{0,2}`
     - `_compute_F_FF__0_2_component <#nachos.core.shaking.Shaker._compute_F_FF__0_2_component>`_
   * -
     - :math:`[\mu^3]^{1,0}`
     - `_compute_F_F_F__1_0_component <#nachos.core.shaking.Shaker._compute_F_F_F__1_0_component>`_
   * -
     - :math:`[\mu^3]^{0,1}`
     - `_compute_F_F_F__0_1_component <#nachos.core.shaking.Shaker._compute_F_F_F__0_1_component>`_
   * - Second hyperpolarizability
     - :math:`[\alpha^2]^{0,0}`
     - `_compute_FF_FF__0_0_component <#nachos.core.shaking.Shaker._compute_FF_FF__0_0_component>`_
   * -
     - :math:`[\alpha^2]^{1,1}`
     - `_compute_FF_FF__1_1_component <#nachos.core.shaking.Shaker._compute_FF_FF__1_1_component>`_
   * -
     - :math:`[\alpha^2]^{2,0}`
     - `_compute_FF_FF__2_0_component <#nachos.core.shaking.Shaker._compute_FF_FF__2_0_component>`_
   * -
     - :math:`[\alpha^2]^{0,2}`
     - `_compute_FF_FF__0_2_component <#nachos.core.shaking.Shaker._compute_FF_FF__0_2_component>`_
   * -
     - :math:`[\mu\beta]^{0,0}`
     - `_compute_F_FFF__0_0_component <#nachos.core.shaking.Shaker._compute_F_FFF__0_0_component>`_
   * -
     - :math:`[\mu\beta]^{1,1}`
     - `_compute_F_FFF__1_1_component <#nachos.core.shaking.Shaker._compute_F_FFF__1_1_component>`_
   * -
     - :math:`[\mu\beta]^{2,0}`
     - `_compute_F_FFF__2_0_component <#nachos.core.shaking.Shaker._compute_F_FFF__2_0_component>`_
   * -
     - :math:`[\mu\beta]^{0,2}`
     - `_compute_F_FFF__0_2_component <#nachos.core.shaking.Shaker._compute_F_FFF__0_2_component>`_
   * -
     - :math:`[\mu^4]^{1,1}`
     - `_compute_F_F_F_F__1_1_component <#nachos.core.shaking.Shaker._compute_F_F_F_F__1_1_component>`_
   * -
     - :math:`[\mu^4]^{2,0}`
     - `_compute_F_F_F_F__2_0_component <#nachos.core.shaking.Shaker._compute_F_F_F_F__2_0_component>`_
   * -
     - :math:`[\mu^4]^{0,2}`
     - `_compute_F_F_F_F__0_2_component <#nachos.core.shaking.Shaker._compute_F_F_F_F__0_2_component>`_
   * -
     - :math:`[\mu^2\alpha]^{1,0}`
     - `_compute_F_F_FF__1_0_component <#nachos.core.shaking.Shaker._compute_F_F_FF__1_0_component>`_
   * -
     - :math:`[\mu^2\alpha]^{0,1}`
     - `_compute_F_F_FF__0_1_component <#nachos.core.shaking.Shaker._compute_F_F_FF__0_1_component>`_


The formulas are detailed :download:`in this document <./formulas_tex/contribs.pdf>`.


.. automodule:: nachos.core.shaking
    :members:
    :private-members:


Analyzing (``nachos.core.analyzing``)
-------------------------------------

.. automodule:: nachos.core.analyzing
    :members: