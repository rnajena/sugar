# (C) 2024, Tom Eulenfeld, MIT license

from importlib.resources import files
import json
from os.path import exists


def _included_color_schemes():
    """
    Return a list of included color schemes
    """
    return sorted(p.stem for p in
                  files('sugar.imaging.color_schemes').joinpath('').glob('*.json'))


def get_color_scheme(name):
    """
    Return color scheme

    :param name: Name of the color scheme,
        the following color schemes are supported out of the box: ``{}``.
        These color schemes originate from Biotite (Gecos), Jaview and ClustalX,
        see `here`_ for a description.
        Alternatively, a JSON file, as generated by the `Gecos`_ library can be specified.
        A ``colors`` entry with a letter->color mapping in the JSON file is sufficient.
    :return: Dictionary with letter->color mapping

    .. _here: https://www.biotite-python.org/latest/examples/gallery/sequence/misc/color_schemes_protein.html#sphx-glr-examples-gallery-sequence-misc-color-schemes-protein-py
    .. _Gecos: https://gecos.biotite-python.org
    """
    if exists(name):
        fname = name
    else:
        fname = str(files('sugar.imaging.color_schemes').joinpath(name.lower() + '.json'))
        if not exists(fname):
            raise ValueError(f'"{name}" not a valid color scheme or file')
    with open(fname) as f:
        cs = json.load(f)
    return cs['colors']


if hasattr(get_color_scheme, '__doc__'):
    get_color_scheme.__doc__ = get_color_scheme.__doc__.format(', '.join(_included_color_schemes()))
