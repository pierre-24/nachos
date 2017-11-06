import os


class BadCooking(Exception):
    pass


class Cooker:
    """Cooker class to retrieve the information out of the calculation results

    :param recipe: a recipe
    :type recipe: nachos.core.files.Recipe"""

    def __init__(self, recipe, directory='.'):
        self.recipe = recipe

        if not os.path.isdir(directory):
            raise BadCooking('{} is not a directory'.format(directory))

        self.directory = directory

    def cook(self):
        """Do the thing"""

        pass
