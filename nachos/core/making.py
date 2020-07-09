import os

from prompt_toolkit import document as pt_document, prompt as prompt_pt
from prompt_toolkit.completion import Completer, Completion
from prompt_toolkit.contrib.completers import WordCompleter, PathCompleter
from prompt_toolkit.validation import Validator, ValidationError
from prompt_toolkit.contrib.validators.base import SentenceValidator

from qcip_tools.chemistry_files import helpers, PropertyNotPresent, PropertyNotDefined, gaussian

from nachos.core import files, CONFIG


class SentenceValidatorWithDefault(SentenceValidator):

    def __init__(self, default=None, **kwargs):
        super().__init__(**kwargs)

        self.default_sentence = default

        if default is not None and default not in self.sentences:
            raise ValueError('default value {} not in list'.format(default))

    def validate(self, document):

        if document.text == '' and self.default_sentence is not None:
            document.text = self.default_sentence
        else:
            super().validate(document)


class ChemistryFileValidator(Validator):
    def __init__(self):
        super().__init__()

        self.validated_geometry = None

    def validate(self, document):
        path = document.text

        if not os.path.exists(path):
            raise ValidationError(message='{} is not a valid file'.format(path))

        try:
            f = helpers.open_chemistry_file(open(path))
        except helpers.ProbablyNotAChemistryFile:
            raise ValidationError(message='unable to open {} as a chemistry file'.format(path))

        try:
            self.validated_geometry = f.property('molecule')
        except (PropertyNotPresent, PropertyNotDefined):
            raise ValidationError(message='unable to find a geometry in {}'.format(path))


class TypeValidator(Validator):
    type = 'None'

    def __init__(self, default=None):
        super().__init__()
        self.converted_value = None
        self.default = default

    def validate(self, document):
        if document.text != '':
            try:
                self.converted_value = self.convert(document.text)
            except ValueError:
                raise ValidationError(message='{} is not a {}'.format(document.text, self.type))
        elif self.default is not None:
            self.converted_value = self.default
        else:
            raise ValidationError(message='empty value is not allowed')

    def convert(self, text):
        raise NotImplementedError()


class TypeIntValidator(TypeValidator):
    type = 'int'

    def convert(self, text):
        return int(text)


class TypeFloatValidator(TypeValidator):
    type = 'float'

    def convert(self, text):
        return float(text)


class BasisCompleter(Completer):

    def __init__(self, bases):
        super().__init__()
        self.bases = bases

    def get_completions(self, document, complete_event):
        text = document.text[:document.cursor_position]
        prev_comma = text.rfind(';')
        current_base = text[prev_comma + 1:]

        if text[-1] != ':':
            for b in self.bases:
                if b.startswith(current_base):
                    yield Completion(b, (prev_comma - len(text)) + 1 if prev_comma > 0 else (-len(text)))

        else:
            current_base = current_base[:-1]
            if current_base in self.bases:
                for i in range(1, 6):
                    yield Completion(str(i), 0)


class BasisValidator(Validator):

    def __init__(self, bases, max_diff, method):
        super().__init__()

        self.bases = bases
        self.max_diff = max_diff
        self.validated_bases = []
        self.method = method

    def validate(self, document):
        text = document.text
        l_text = len(text)

        if l_text == 0:
            raise ValidationError(message='differentiation is empty')

        self.validated_bases = []
        cursor_location = 0
        already_in = []

        for x in text.split(';'):
            info = x.split(':')
            cursor_location += len(x)

            if len(info) != 2:
                raise ValidationError(
                    message='{} is not a correct input (should be "basis:level")'.format(x),
                    cursor_position=cursor_location)
            else:
                b_repr = info[0]
                try:
                    level = int(info[1])
                except ValueError:
                    raise ValidationError(
                        message='level of differentiation {} (associated with {}) is not a number'.format(
                            info[1], info[0]),
                        cursor_position=cursor_location)

                if level < 1:
                    raise ValidationError(
                        message='level of differentiation for {} is too small (should be >0)'.format(b_repr),
                        cursor_position=cursor_location)

            if b_repr in already_in:
                raise ValidationError(
                    message='{} is defined multiple times'.format(b_repr),
                    cursor_position=cursor_location)

            if b_repr not in self.bases:
                raise ValidationError(
                    message='{} is not allowed for this flavor'.format(b_repr),
                    cursor_position=cursor_location)

            num_G = b_repr.count('G')
            num_F = b_repr.count('F') + b_repr.count('D')

            if num_G > self.max_diff['G'] or num_F > self.max_diff['F']:
                raise ValidationError(
                    message='"{}" for method {} is not possible'.format(b_repr, self.method),
                    cursor_position=cursor_location)

            self.validated_bases.append((b_repr, level))
            already_in.append(b_repr)

            cursor_location += 1  # after next column

        if len(self.validated_bases) == 0:
            raise ValidationError(message='differentiation is empty')

    def differentiation(self):
        diff = {}
        for basis, level in self.validated_bases:
            if level not in diff:
                diff[level] = []

            diff[level].append(basis)

        return diff


class FrequenciesValidator(Validator):

    def __init__(self):
        super().__init__()
        self.validated_frequencies = []

    def validate(self, document):
        text = document.text
        self.validated_frequencies = []
        cursor_location = 0

        allowed_units = ['nm', 'ev', 'cm-1']

        for x in text.split(';'):
            cursor_location += len(x)

            value = None

            try:
                value = float(x)
            except ValueError:
                for unit in allowed_units:
                    if x[-len(unit):] == unit:
                        try:
                            float(x[:-len(unit)])
                            value = x
                        except ValueError:
                            pass
                        break

            if value is None:
                raise ValidationError(
                    message='{} is not an allowed frequency'.format(x), cursor_position=cursor_location)

            self.validated_frequencies.append(value)
            cursor_location += 1


class GenBasisValidator(Validator):

    def __init__(self, geometry):
        super().__init__()
        self.geometry = geometry

    def validate(self, document):
        text = document.text
        gbs = gaussian.BasisSet()

        if not os.path.exists(text):
            raise ValidationError(message='cannot open {}'.format(text))

        try:
            gbs.read(open(text))
        except gaussian.BasisSetFormatError as e:
            raise ValidationError(message='error while opening basis set in {}: {}'.format(text, str(e)))

        for i in self.geometry.symbols_contained:
            if i not in gbs.basis_set:
                raise ValidationError(message='No basis set for {} defined'.format(i))


class ExtraFlavorCompleter(Completer):
    def __init__(self, allowed_keywords):
        super().__init__()
        self.allowed_keywords = allowed_keywords

    def get_completions(self, document, complete_event):
        text = document.text[:document.cursor_position]
        prev_comma = text.rfind(';')
        current_base = text[prev_comma + 1:]

        if text[-1] != '=':
            for b in self.allowed_keywords:
                if b.startswith(current_base):
                    yield Completion(b, (prev_comma - len(text)) + 1 if prev_comma > 0 else (-len(text)))


class ExtraFlavorValidator(Validator):

    def __init__(self, allowed_keywords):
        super().__init__()
        self.allowed_keywords = allowed_keywords
        self.validated_keywords = {}

    def validate(self, document):
        text = document.text

        self.validated_keywords = {}
        cursor_location = 0

        if len(text) == 0:
            return

        for x in text.split(';'):
            info = x.split('=')
            cursor_location += len(x)

            if len(info) < 2:
                raise ValidationError(
                    message='{} is not a correct input (should be "var=value")'.format(x),
                    cursor_position=cursor_location)
            elif len(info) > 2:
                info[1] = '='.join(info[1:])
                info = info[:2]

            keyword, value = info
            if keyword not in self.allowed_keywords:
                raise ValidationError(
                    message='{} is not an allowed keyword'.format(keyword),
                    cursor_position=cursor_location)

            self.validated_keywords[keyword] = value
            cursor_location += 1


class BadMaking(Exception):
    pass


class MakingAction:
    def __init__(self):
        pass

    def use(self, variable, value, validator=None):
        raise NotImplementedError('use')


class SetRecipeAction(MakingAction):
    def __init__(self, recipe):
        super().__init__()
        self.recipe = recipe

    def use(self, variable, value, validator=None):
        kw = {variable: value}
        self.recipe.recipe.update(**kw)
        self.recipe._update(kw)


class SetRecipeGeometryAction(SetRecipeAction):
    def use(self, variable, value, validator=None):
        if validator is None:
            raise ValueError('validator should not be empty')

        self.recipe['geometry'] = value
        self.recipe.geometry = validator.validated_geometry
        self.recipe.dof = 3 * len(validator.validated_geometry)


class SetRecipeBasesAction(SetRecipeAction):
    def use(self, variable, value, validator=None):
        if validator is None:
            raise ValueError('validator should not be empty')

        super().use(variable, validator.differentiation())


class SetRecipeWithConversionAction(SetRecipeAction):
    def use(self, variable, value, validator=None):
        if validator is None:
            raise ValueError('validator should not be empty')
        if not isinstance(validator, TypeValidator):
            raise ValueError('validator must be TypeValidator')

        super().use(variable, validator.converted_value)


class SetRecipeExtraFlavorAction(SetRecipeAction):
    def use(self, variable, value, validator=None):
        kw = {variable: value}
        self.recipe.recipe['flavor_extra'].update(**kw)


class SetRecipeAllExtraFlavorAction(SetRecipeAction):
    def use(self, variable, value, validator=None):
        if validator is None:
            raise ValueError('validator should not be empty')
        self.recipe.recipe['flavor_extra'].update(**validator.validated_keywords)


class SetRecipeFrequenciesAction(SetRecipeAction):
    def use(self, variable, value, validator=None):
        if validator is None:
            raise ValueError('validator should not be empty')

        super().use(variable, validator.validated_frequencies)


class Maker:
    """Make a recipe out of prompting question or arguments"""

    def __init__(self, use_fallback_prompt=False, raise_when_arg_wrong=False):
        self.use_fallback_prompt = use_fallback_prompt
        self.raise_when_arg_wrong = raise_when_arg_wrong

    def make(self, args):
        """make the recipe

        :param args: input arguments
        :type args: dict
        :rtype: nachos.core.files.Recipe
        """

        recipe = files.Recipe()

        self._make_var(
            args,
            'flavor',
            'What flavor for you, today?',
            SetRecipeAction(recipe),
            completer=WordCompleter(words=CONFIG.keys()),
            validator=SentenceValidatorWithDefault(sentences=CONFIG.keys())
        )

        config = CONFIG[recipe['flavor']]
        self._make_var(
            args,
            'type',
            'What type of differentiation?',
            SetRecipeAction(recipe),
            completer=WordCompleter(words=config['types']),
            validator=SentenceValidatorWithDefault(sentences=config['types']))

        methods = [a[0] for a in config['methods']]
        self._make_var(
            args,
            'method',
            'With which method?',
            SetRecipeAction(recipe),
            completer=WordCompleter(words=methods),
            validator=SentenceValidatorWithDefault(sentences=methods))

        if recipe['method'] == 'DFT':
            self._make_var(
                args,
                'XC',
                'Which XC functionnal?',
                SetRecipeExtraFlavorAction(recipe))

        if recipe['method'] == 'CC':
            self._make_var(
                args,
                'CC',
                'Which Coupled Cluster level?',
                SetRecipeExtraFlavorAction(recipe))

        max_diff = next(a[1] for a in config['methods'] if a[0] == recipe['method'])

        self._make_var(
            args,
            'geometry',
            'Where is the geometry?',
            SetRecipeGeometryAction(recipe),
            completer=PathCompleter(),
            validator=ChemistryFileValidator())

        self._make_var(
            args,
            'basis_set',
            'With which basis set?',
            SetRecipeAction(recipe))

        if recipe['basis_set'] == 'gen':
            self._make_var(
                args,
                'gen_basis',
                'Where is the gen basis set?',
                SetRecipeExtraFlavorAction(recipe),
                completer=PathCompleter(),
                validator=GenBasisValidator(recipe.geometry))

        self._make_var(
            args,
            'differentiation',
            'What to differentiate?',
            SetRecipeBasesAction(recipe),
            completer=BasisCompleter(config['bases']),
            validator=BasisValidator(config['bases'], max_diff, recipe['method']))

        need_frequency = False
        for level in recipe['differentiation']:
            for basis in recipe['differentiation'][level]:
                if 'D' in basis:
                    need_frequency = True
                    break
                if need_frequency:
                    break

        if need_frequency:
            self._make_var(
                args,
                'frequencies',
                'Dynamic frequencies?',
                SetRecipeFrequenciesAction(recipe),
                validator=FrequenciesValidator())

        self._make_var(
            args,
            'name',
            'Name of the files?',
            SetRecipeAction(recipe),
            default=files.DEFAULT_RECIPE['name'])

        self._make_var(
            args,
            'min_field',
            'Minimum field (F0)?',
            SetRecipeWithConversionAction(recipe),
            validator=TypeFloatValidator(default=files.DEFAULT_RECIPE['min_field']),
            default=files.DEFAULT_RECIPE['min_field'])

        self._make_var(
            args,
            'ratio',
            'Ratio (a)?',
            SetRecipeWithConversionAction(recipe),
            validator=TypeFloatValidator(default=files.DEFAULT_RECIPE['ratio']),
            default=files.DEFAULT_RECIPE['ratio'])

        self._make_var(
            args,
            'k_max',
            'Maximum k?',
            SetRecipeWithConversionAction(recipe),
            validator=TypeIntValidator(default=files.DEFAULT_RECIPE['k_max']),
            default=files.DEFAULT_RECIPE['k_max'])

        self._make_var(
            args,
            'flavor_extra',
            'Update flavor extra ? (left blank for no change)',
            SetRecipeAllExtraFlavorAction(recipe),
            validator=ExtraFlavorValidator(config['default_for_extra_fields']),
            completer=ExtraFlavorCompleter(config['default_for_extra_fields']))

        return recipe

    def _make_var(self, args, variable, message, action, validator=None, completer=None, default=None):
        """Make the action for a give variable

        :param args: input args
        :type args: dict
        :param variable: variable
        :type variable: str
        :param message: prompt message:
        :type message: str
        :param validator: the validator, if any
        :type validator: prompt_toolkit.validation.Validator
        :param completer: completer, if any (relevant for prompt_toolkit)
        :type completer: prompt_toolkit.completion.Completer
        :param action: the action to perform
        :type action: MakingAction
        :param default: default value, if any
        :type default: str
        """

        value = None

        if hasattr(args, variable):
            value = getattr(args, variable, default)

            if validator and value:
                try:
                    validator.validate(pt_document.Document(text=value))
                except ValidationError as e:
                    if self.raise_when_arg_wrong:
                        raise BadMaking('error while validating {}: {}'.format(variable, e.message))
                    else:
                        print('error while validating {}: {}'.format(variable, e.message))
                    value = None

        if value is None:
            value = self.prompt(message, validator, completer, default)

        action.use(variable, value, validator)

    def prompt(self, message, validator, completer, default=None):
        if self.use_fallback_prompt:
            return Maker._prompt_fallback(message, validator, completer, default)
        else:
            return Maker._prompt_toolkit(message, validator, completer, default)

    @staticmethod
    def _prompt_fallback(message, validator, completer, default=None):
        while True:
            r = input('{}{} '.format(message, ' [{}]'.format(default) if default else ''))

            if r == '' and default:
                r = default

            if validator:
                try:
                    validator.validate(pt_document.Document(text=str(r)))
                except ValidationError as e:
                    print('error: {}'.format(e.message))
                    continue
                break
            else:
                break
        return r

    @staticmethod
    def _prompt_toolkit(message, validator, completer, default=None):

        v = prompt_pt(
            message='{}{} '.format(message, ' [{}]'.format(default) if default is not None else ''),
            validator=validator,
            completer=completer
        )

        if default is not None and validator is None and v == '':
            v = default

        return v
