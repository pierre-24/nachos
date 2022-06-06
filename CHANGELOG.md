# Changelog

## Current version

+ Add inversion in analyze
+ Group vibrational contributions
+ Add more information in `nachos_cook` and catch correct energy with Gaussian
+ Add Q-Chem, but only for [CCMAN2](http://www.q-chem.com/qchem-website/manual/qchem51_manual/sect-ccmeth.html) energies
+ Output version and parameters (as requested by Benoît)
+ Add uncertainties to the output.
+ Switch to `pip-tools`
+ Go to Github, and use `bump2version` and GH actions.
+ Add SCS-MP2 and ability to extract data from Gaussian LOG file.
+ Skip the auto procedure for choosing best value in Romberg triangle.
+ Add `nachos_peek` to look into files.

## Version 0.3

+ Upgrade to latest version of qcip_tools (0.5.3.2)
+ Multiple input directories for cooking
+ Add order II pv contributions to second hyperpolarizability (`[α²]¹¹`, `[µβ]¹¹`, `[µ⁴]¹¹`, `[α²]²⁰`, `[α²]⁰²`, `[µβ]²⁰`, `[µβ]⁰²`, `[µ⁴]²⁰`, `[µ⁴]⁰²`, finish #4)
+ Add a frequency limiter to `nachos_analyze`

## Version 0.2

+ In arguments, lists accepts empty elements
+ Separate dalton *dal* inputs
+ Allow to perform HF and DFT calculation with dalton
+ Allow to contribute pv contributions if only the static properties are available
+ Add order 0 and I pv contributions to second hyperpolarizability (`[α²]⁰⁰`, `[µβ]⁰⁰`, `[µ²α]¹⁰`, `[µ²α]⁰¹`)
+ Upgrade to latest version of qcip_tools (0.5.3).

## Version 0.1.1

+ Bugfix: a too strong assumption on the symmetry of the tensor leads to incorrect values for certain components.
+ Add a changelog.

## Version 0.1

First version. All program created and working.

ZPVA available for any property. Pure vibrational contributions available for the polarizability and the first hyperpolarizability.