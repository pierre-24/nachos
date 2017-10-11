# Contributing to NACHOS

Please note that the code is not actually developed on the git server of the University of Namur (which only contains the releases) but on a personal protected git server (with CI activated and properly configured).
Feel free to ask access if needed.

First, you need to [install](README.md#installation) the repository.

Maybe you also want to take a look into [qcip_tools](https://gitlab.unamur.be/pierre.beaujean/qcip_tools), the library behind the scripts.

## Design rules

+ The code is written in Python 3, and follows the (in)famous [PEP-8 ](http://legacy.python.org/dev/peps/pep-0008/). You can check it by running `make lint`, which launch the `flake` utility.
+ Codes and comments are written in english.
+ If you add a new script, please edit `provided_scripts` in `nachos/__init__py` to add yours, so that it is added to the installed scripts.
+ The code is tested. 
  You can launch the test series by using `make test`. 
  Every script should be provided with at least one unit test. 
  You may need test files to do so, but try to make them small (say, don't use d-aug-cc-pVDZ while STO-3G could do the job).
  But there is no need to check that it works for every kind of file, since it relies on the property on the `qcip_tools` side!
  Rather than that, use a file that works, one that don't, and test program options.


## Workflow

Adapted from the (in)famous [Git flow](http://nvie.com/posts/a-successful-git-branching-model/).

+ Development is mad in `dev` branch, while `master` contains the production version (and is protected from editing).
+ Functionalities are added through merge request (MR) in the `dev` branch. Do not work in `dev` directly, but create a new branch (`git checkout -b my_branch origin/dev`).
+ Theses merge requests should be unitary, and include unit test(s) and documentation if needed. The test suite must succeed for the merge request to be accepted.
+ At some (random) points, `dev` will be merged by the maintainer into `master` to create a new version, with a tag of the form `release-vXX`.

## Licence

This code belong to me, [Pierre Beaujean](pierre.beaujean@unamur.be), and to the [University of Namur](https://www.unamur.be) since it is developed and used in the frame of my PhD thesis.
Of course, if you contribute, you will be added to this list ;)
