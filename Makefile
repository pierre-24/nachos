all: help

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  install-dependencies        to install python dependencies through pip"
	@echo "  install-dependencies-dev    to install python dependencies (for dev) through pip"
	@echo "  lint                        to lint backend code (flake8)"
	@echo "  test                        to run test suite"
	@echo "  help                        to get this help"

install-dependencies:
	pip3 install --upgrade -r requirements.txt

install-dependencies-dev: install-dependencies
	pip3 install --upgrade -r requirements-dev.txt

lint:
	flake8 nachos tests --max-line-length=120 --ignore=N802

test:
	python setup.py test