=================
How to use nachos
=================

Concepts
--------

+  The numerical differentiation is used to obtain (static) derivatives to different quantities.
   It exploits the following identity, given :math:`f(x)` a given function expanded in Taylor series,

   .. math::

      \left.\frac{\partial f(x)}{\partial x}\right|_{x=0} = \lim_{h_0\rightarrow 0} \frac{f(h_0)-f(0)}{h_0},

   where :math:`h_0` is the minimal field (``min_field`` parameter).

+  The code use the Romberg procedure to remove contamination from higher orders (see `this publication <dx.doi.org/10.1002/qua.24685>`_ for more details).
   The derivative is computed for different values of :math:`f=a^k\,f_0`, with :math:`k<k_{max}` the field amplitude (which lead to the ``k_max`` parameter), and :math:`a` is the common ratio (``ratio`` parameter).
   The procedure goes as follow:

   .. math::

      \begin{align}
      &H_{k,0} = \frac{f(a^k\,h_0)-f(0)}{a^k\,h_0},\\
      &H_{k,m+1} = \frac{a^{2m}\,H_{k,m}-H_{k+1,m}}{a^{2m}-1},
      \end{align}

   where :math:`m` is the number of iterations (or refinement steps).
   This leads to a so-called *Romberg triangle*, from which the value of the derivative is extracted.

+  The different quantities are written **as derivatives with respect to the energy**.
   For example, the geometric Hessian, second order of the energy with respect to geometrical derivatives, is written ``GG``, ``G`` meaning *geometrical derivative with respect to cartesian coordinates*.
   Accordingly,

   + ``F`` means *derivatives with respect to static electric field* ;
   + ``D`` means *derivatives with respect to dynamic electric field* (with a given frequency), and ``d`` means the same, but with inverse frequency (:math:`-\omega`) ;
   + ``N`` means *derivatives with respect to normal coordinates*.

   Therefore, the static hyperpolarizability, :math:`\beta(0;0,0)`, is written ``FFF``, while the dynamic hyperpolarizability depends on the process involved: ``FDF`` for EOP  [:math:`\beta(-\omega;\omega,0)`] and ``FDD`` for SHG [:math:`\beta(-2\omega;\omega,\omega)`].
   See the list `below <#list-of-the-derivatives>`_.

   Geometrical derivatives of an electrical derivative are written with the geometrical derivatives **first**.
   For example, the first order geometrical derivative (with respect to normal mode) of the static polarizability, :math:`\frac{\partial \alpha}{\partial Q}`, is written ``NFF``, the second order one ``NNFF``.

   Note that the number of ``G`` and ``N`` thus correspond to the level of geometric differentiation and the number of ``F``, ``D`` and ``d`` to the level of electrical differentiation.

+  Nachos is abble to perform differentiation with respect to static electric field (``F``) and cartesian coordinate (``G``).
   Given the cartesian hessian, ``nachos_bake`` (see below) perform a vibrational analysis and is able to project ``G`` derivatives over normal mode, giving the corresponding ``N`` ones.


General workflow
----------------

Here is the schematic of the workflow with the nachos package:


.. figure:: ./images/workflow.png
   :align: center

   Flowchart for the different parts of the nachos package. Arrows indicate whether a part is an input (arrows going in) of a program (rectangle) or an output (arrow going out).

In short,

1. ``nachos_make`` create a *recipe* (``nachos_recipe.yml``, but you can change that), which is the file that explains what to do and how to do it ;
2. ``nachos_prepare`` uses the recipe to generate the different input files for the quantum chemistry program of your choice (currently Gaussian and Dalton) ;
3. The quantum chemistry program process the different input files and generate output files ;
4. ``nachos_cook`` carry out all the information that it can get from output files (FCHK files for gaussian, TAR archives and OUT files for Dalton) and store them in a *data file* (``nachos_data.h5``, but you can change that) ;
5. ``nachos_bake`` perform the requested numerical differentiation(s) out of the data from the *data file*, and store them in a *final file* (``molecule_nd.h5``, but you can change that) ;
6. (*optional*) ``nachos_analyze`` allow to quickly get a given property for each quantity stored in the *final file* (e.g. a tensor component, an average, ...) ;
7. (*optional*) if possible, ``nachos_shake`` will add the different vibrational contribution to the electrical derivatives of the energy.

Therefore, it should looks like:

.. sourcecode:: bash

    # create subdirectory:
    mkdir new_directory
    cd new_directory

    # create recipe and inputs:
    nachos_make && nachos_prepare

    # ... run all the inputs ...

    # carry out and perform numerical differentiation:
    nachos_cook && nachos_bake # or nachos_bake -P

    # ... eventually, add vibrational contributions:
    nachos_shake

    # ... eventually, post-analyze:
    nachos_analyze -p "xxxx" > nachos.log


See below for more details on every command.

.. autoprogram:: nachos.make:get_arguments_parser()
    :prog: nachos_make

.. note::

    + It is easier to place the geometry file (and eventual basis set and other extra files) in the **same** directory as the recipe.
    + For some terminal, it is not possible to use the extended prompt toolkit, use ``-N`` to get an alternative.
    + Default behavior is if there is an error in the input argument, the corresponding question is asked again.
      If you just want the program to fail (because you are using it in a script), use the ``-S`` option.
    + ``F`` differentiation is **only possible** with gaussian.

The program prompts for different information in order to create a *recipe file*, if not given in command line, and generate a recipe in output (``-o`` option, default is ``nachos_recipe.yml``).

.. list-table::
   :header-rows: 1
   :widths: 20 35 35 10

   * - Option
     - Question
     - Possible inputs
     - Note
   * - ``--flavor``
     - "What flavor for you, today?"
     - ``gaussian`` | ``dalton``
     -
   * - ``--type``
     - "What type of differentiation?"
     - ``F`` | ``G``
     -
   * - ``--method``
     - "With which method?"
     - :ref:`see below <nachos_make_note_1>`
     -
   * - ``--XC``
     - "Which XC functionnal?"
     - *XC functional*
     - Only if ``DFT``
   * - ``--geometry``
     - "Where is the geometry? "
     - *path to a .com/.xyz/.fchk/.mol* file
     -
   * - ``--basis-set``
     - "With which basis set?"
     - *valid basis set* | ``gen``
     -
   * - ``--gen-basis``
     - "Where is the gen basis set?"
     - *path to a gbs file*
     - Only if ``gaussian`` and ``gen``
   * - ``--differentiation``
     - "What to differentiate?"
     - :ref:`see below <nachos_make_note_2>`
     -
   * - ``--frequencies``
     - "Dynamic frequencies?"
     - :ref:`see below <nachos_make_note_3>`
     - Only if dynamic quantities requested
   * - ``--name``
     - "Name of the files?"
     - *any string*
     - Avoid spaces and special characters!
   * - ``--min-field``
     - "Minimum field (F0)?"
     - *floating number*
     -
   * - ``--ratio``
     - "Ratio (a)?"
     - *floating number*
     -
   * - ``--k-max``
     - "Maximum k?"
     - *floating number*
     -
   * - ``--flavor-extra``
     - "Update flavor extra ?"
     - :ref:`see below <nachos_make_note_4>`
     - Blank input use default values

When everything is done, you end up with a ``.yml`` file that contains all the information you input.
For example, this is an input to compute vibrational contribution to the polariability:

.. code-block:: yaml

    # flavor
    flavor: gaussian
    method: HF
    basis_set: gen
    geometry: water.xyz
    flavor_extra:
      convergence: 11
      cphf_convergence: 10
      gen_basis: sto-3g.gbs
      memory: 3Gb
      procs: 4
    # differentiation (the label is the number of time
    # you want to differentiate each item of the list)
    differentiation:
      2:
        - F
        - FF
        - FD
      1:
        - GG
    type: G
    min_field: 0.01
    ratio: 2
    k_max: 3
    frequencies:
      - 1064nm
      - 694.3nm
    # others:
    name: water_test

Obviously, nothing prevents you from writing your own *recipe file* from scratch. Actually, you just need to define

    + ``flavor`` ;
    + ``type`` ;
    + ``method`` ;
    + ``basis_set`` ;
    + ``geometry`` ;
    + ``differentiation`` ;

Since there is default values for the rest.

-------

.. _nachos_make_note_1:

For ``--method``: the value of this argument depends on the *flavor* you chose.
This also determine the maximum derivative available at this level i.e. what you can request in ``--differentiation`` (:ref:`see below <nachos_make_note_2>`).

+ For ``gaussian`` (chosen according to the `force page <http://gaussian.com/force/>`_, the `freq page <http://gaussian.com/freq/>`_ and the `polar page <http://gaussian.com/polar/>`_):

  .. list-table::
       :header-rows: 1
       :widths: 30 20 20 30

       * - Method
         - Maximum level of electrical differentiation
         - Maximum level of geometrical differentiation
         - Available
       * - ``HF``
         - 3
         - 2
         - ``energy``, ``G``, ``GG``, ``F``, ``FF``, ``FD``, ``FDF``, ``FDD``
       * - ``DFT``
         - 3
         - 2
         - ``energy``, ``G``, ``GG``, ``F``, ``FF``, ``FD``, ``FDF``, ``FDD``
       * - ``MP2``
         - 2
         - 2
         - ``energy``, ``G``, ``GG``, ``F``, ``FF``
       * - ``MP3``, ``MP4``, ``MP4D``, ``MP4DQ``, ``MP4SDQ``
         - 1
         - 1
         - ``energy``, ``G``, ``F``
       * - ``CCSD``
         - 1
         - 1
         - ``energy``, ``G``, ``F``
       * - ``CCSD(T)``
         - 0
         - 0
         - ``energy``

  Some method are not available, but may be added in the future if needed (CI methods, for example).

+ For ``dalton`` you can request ``CCS``, ``CC2``, ``CCSD`` and ``CC3``, for which you can request derivatives up to second hyperpolarizability, and the gradient.

-------

.. _nachos_make_note_2:

For ``--differentiation``: this is where you request what you want to differentiate, and up to which level, with a semicolon separated list.
Each member of the list should be of the form ``what:how many``, where ``what`` is a derivative (`see the appendix <#list-of-derivatives>`_) and ``how much`` is how many times you want to differentiate this quantity.

For example,

+ If you want to do an electric field differentiation (``F``) to obtain the static first hyperpolarizability (``FFF``) from the energy, input should be ``energy:3``, because you want to differentiate energy 3 times.
  To get the same property from the dipole moment and the static polarizability, the input is ``F:2;FF:1``.
+ If you want to get the vibrational contribution to a given property (say, the polarizability), you need to select ``G`` for the type of differentiation, then you need at least second order derivative of the dipole moment polariability with respect to that (the first one is automatically computed if the second is), and the cubic force field, so an input could look like ``FF:2;F:2;GG:1`` (and eventually ``FD:2``).

:ref:`See above <nachos_make_note_1>` for the list of quantities that you can differentiate depending on the *flavor* and the method.

-------

.. _nachos_make_note_3:

For ``--frequencies``: This is only relevant if you requested the differentation of a quantity that is dynamic.
The input is a list of semicolon separated frequencies, and is quite liberal, since a valid example could be ``1064nm;0.04:1000cm-1;0.1eV`` (it accepts ``eV``, ``cm-1``,  ``nm`` and nothing, which means atomic units).
The values are converted in atomic unit in ``nachos_prepare`` (see below).

-------

.. _nachos_make_note_4:

For ``--flavor-extra``: this option actually controls the generation of input files and that is it (for example, that is where you request the amount of memory and processors for gaussian).
The options depends on the *flavor*, and are given in a semicolon separated list (for example ``procs=4;memory=3Gb;extra_keywords=srcf=(iefpcm,solvent=water)`` for ``gaussian``).
Note that you don't have to redefine every variable, since they have a default value which is correct for most cases.

+ For ``gaussian``, the options are

  .. list-table::
       :header-rows: 1
       :widths: 20 20 60

       * - Option
         - Default value
         - Note
       * - ``memory``
         - ``1Gb``
         - Value of ``%mem``
       * - ``procs``
         - ``1``
         - Value of ``%nprocshared``
       * - ``convergence``
         - ``11``
         - SCF convergence criterion
       * - ``cphf_convergence``
         - ``10``
         - CPHF convergence criterion
       * - ``cc_convergence``
         - ``11``
         - CC convergence criterion
       * - ``max_cycle``
         - ``600``
         - Maximum number of SCF and CC cycles
       * - ``extra_keywords``
         -
         - Any extra input (for example, the solvent, ...)
       * - ``extra_section``
         -
         - Path to a file where extra section of the input files are given (for example, solvent definition, ...)
       * - ``vshift``
         - ``1000``
         - Apply a *vshift* (helps for the electric field differentiation)

  Note that the value of ``extra_section`` is not tested here.
  Also, ``XC`` and ``gen_basis`` are available, but that would increase their previous values.

+ For ``dalton``, the options are

  .. list-table::
       :header-rows: 1
       :widths: 20 20 60

       * - Option
         - Default value
         - Note
       * - ``max_iteration``
         - ``2500``
         - Maximum number of iteration for the response function computation
       * - ``threshold``
         - ``1e-6``
         - Convergence criterion for the SCF
       * - ``cc_threshold``
         - ``1e-11``
         - Convergence criterion for CC energy and response functions
       * - ``dal_name``
         - ``ND``
         - Prefix for the different ``.dal`` files


.. autoprogram:: nachos.prepare:get_arguments_parser()
    :prog: nachos_prepare


The program will prepare as many input files as needed.
By using ``-d``, you can decide where the input files should be generated, but keep in mind that they should be in the same directory as the recipe for the next step (use ``-c`` if needed).

The ``-V 1`` option allows you to know how much files where generated.

.. note::

    To helps the dalton program, a file called ``inputs_matching.txt`` is created for this *flavor*, where each lines contains the combination of dal and mol file to launch (because there may be different dal files).

    If you use job arrays, you may therefor use a job file that contains the following lines (here with  `slurm <https://slurm.schedmd.com/>`_, but it is the same with other schedulers):

    .. code-block:: bash

      # get the files from the line:
      INPUT_FILES=$(sed -n "${SLURM_ARRAY_TASK_ID}p" inputs_matching.txt)
      # launch dalton:
      dalton $INPUT_FILES

    You need to launch as many calculations as there is lines in this file.


.. autoprogram:: nachos.cook:get_arguments_parser()
    :prog: nachos_cook

The program fetch the different computational results from each files that it can fin (it looks for FCHK files with gaussian, TAR archive and OUT files for dalton), and mix them together in a single *data file*.

The ``-V 1`` option allows you to know which files the program actually discovered and used.


.. warning::

    The program looks for output files **in the same directory as the recipe**, and there is no way to change this behavior.

Appendix
--------

List of the derivatives
***********************

Note that it would be better to respect the order for the different derivatives (``FDF``, not ``FFD``, for example).

.. list-table::
   :header-rows: 1
   :widths: 40 10 50

   * - Derivative
     -
     - Comment
   * - The energy
     - ``energy``
     -
   * - :math:`\mu`
     - ``F``
     - Dipole moment
   * - :math:`\alpha(0;0)`
     - ``FF``
     - Static polarizability
   * - :math:`\alpha(-\omega;\omega)`
     - ``FD``
     - Dynamic polarizability
   * - :math:`\beta(0;0,0)`
     - ``FFF``
     - Static first hyperpolarizability
   * - :math:`\beta(-\omega;\omega,0)`
     - ``FDF``
     - EOP first hyperpolarizability
   * - :math:`\beta(-2\omega;\omega,\omega)`
     - ``FDD``
     - SHG first hyperpolarizability
   * - :math:`\gamma(0;0,0,0)`
     - ``FFFF``
     - Static second hyperpolarizability
   * - :math:`\gamma(-\omega;\omega,0,0)`
     - ``FDFF``
     - Kerr second hyperpolarizability
   * - :math:`\gamma(-2\omega;\omega,\omega,0)`
     - ``FDDF``
     - ESHG second hyperpolarizability
   * - :math:`\gamma(-\omega;\omega,\omega,-\omega)`
     - ``FDDd``
     - DFWM second hyperpolarizability
   * - :math:`\gamma(-3\omega;\omega,\omega,\omega)`
     - ``FDDD``
     - THG second hyperpolarizability
   * - :math:`\frac{\partial V(x)}{\partial x}`
     - ``G``
     - (cartesian) gradient
   * - :math:`\frac{\partial^2 V(x,y)}{\partial x\partial y}`
     - ``GG``
     - (cartesian) hessian