# Changelog

## Current version

+ Add inversion in analyze (!16)
+ Group vibrational contributions (!17)
+ Update to pip 10 (!18)
+ Switch to pipenv (!19) (and upgrade to the latest version of qcip_tools)
+ Add more information in `nachos_cook` and catch correct energy with Gaussian (!22)
+ Add Q-Chem (!23), but only for [CCMAN2](http://www.q-chem.com/qchem-website/manual/qchem51_manual/sect-ccmeth.html) energies
+ Output version and parameters (as requested by Benoît)

## Version 0.3

+ Upgrade to latest version of qcip_tools (0.5.3.2)
+ Multiple input directories for cooking (!13)
+ Add order II pv contributions to second hyperpolarizability (`[α²]¹¹`, `[µβ]¹¹`, `[µ⁴]¹¹`, `[α²]²⁰`, `[α²]⁰²`, `[µβ]²⁰`, `[µβ]⁰²`, `[µ⁴]²⁰`, `[µ⁴]⁰²`, finish #4)
+ Add a frequency limiter to `nachos_analyze` (#9)

## Version 0.2

+ In arguments, lists accepts empty elements (#3)
+ Separate dalton *dal* inputs (#5) 
+ Allow to perform HF and DFT calculation with dalton (#6)
+ Allow to contribute pv contributions if only the static properties are available
+ Add order 0 and I pv contributions to second hyperpolarizability (`[α²]⁰⁰`, `[µβ]⁰⁰`, `[µ²α]¹⁰`, `[µ²α]⁰¹`, partially #4)
+ Upgrade to latest version of qcip_tools (0.5.3).

## Version 0.1.1

+ Bugfix: a too strong assumption on the symmetry of the tensor leads to incorrect values for certain components (#8).
+ Add a changelog (#7) .


## Version 0.1

First version. All program created and working.

ZPVA available for any property. Pure vibrational contributions available for the polarizability and the first hyperpolarizability.