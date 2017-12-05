# NACHOS

Code for numerical differentiation of energy (oriented toward computation of vibrational contributions).
Maintained by [Pierre Beaujean](pierre.beaujean@unamur.be) and created in the frame of my PhD thesis in the [University of Namur](https://www.unamur.be).

Based on [qcip_tools](https://gitlab.unamur.be/pierre.beaujean/qcip_tools) version [0.4.1](https://gitlab.unamur.be/pierre.beaujean/qcip_tools/tree/release-v0.4.1).
## Installation

### To use the scripts

To just use the scripts contained in the package, just use pip:

```bash
# install qcip_tools (base library):
pip3 install --user --upgrade git+ssh://git@git.pierrebeaujean.net/pierre/qcip_tools.git@dev
# install the scripts:
pip3 install --user --upgrade git+ssh://git@git.pierrebeaujean.net/pierre/nachos.git
```

Note that `--user` allow you to install the package without `sudo` (see [here](https://pip.pypa.io/en/stable/user_guide/#user-installs)).
You will probably need to add `$HOME/.local/bin` to `$PATH` for this to work:

```bash
echo 'PATH=$HOME/.local/bin:$PATH' >> ~/.bashrc
```

On the other hand, you can install it in a *virtualenv* (see step 3 below).

### To contribute

To contribute to the package development: 

1. Clone the repository: `git clone git@git.pierrebeaujean.net:pierre/nachos.git`. To update it, use `git pull` from time to time.
2. Go into the directory: `cd nachos`.
3. Use a *virtualenv* if you don't want this to mess up with your main python installation:

    ```bash
    # create the virtualenv with python 3 as the default python:
    virtualenv venv --python=python3
    # activate it, do that every time you want to use it:
    source venv/bin/activate
    ``` 

4. Install python dependencies: `make install-dependencies-dev`.
5. "Install" the package : `python setup.py develop` (make the different commands are available).

See [the contribution page](CONTRIBUTING.md) for more informations.

## Documentation

*None yet* ;)

## Contributing

You can reports bugs and suggestions any time by email or using the bugtracker.

If you want to contribute to the code, see the [contribution page](CONTRIBUTING.md). 
Please note that the code is not actually developed on the git server of the University of Namur (which only contains the releases) but on a personal protected git server (with CI activated and properly configured). 
Feel free to ask access if needed.
