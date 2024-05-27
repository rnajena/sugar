# (C) 2024, Tom Eulenfeld, MIT license

import collections.abc
import copy


class Attr(collections.abc.MutableMapping):
    """
    A class which behaves like a dictionary.

    :param dict data: Dictionary with initial keywords.

    .. rubric:: Basic Usage

    You may use the following syntax to change or access data in this class.

    >>> meta = Meta()
    >>> meta.comment = 'bla'
    >>> meta['another_comment'] = 'yeah'
    >>> print(stats.get('comment'))
    bla
    >>> print(stats['comment'])
    bla
    >>> print(stats.comment)
    bla
    """

    def __init__(self, *args, **kwargs):
        """
        A Meta object can be initialized in two ways. It can be given an
        existing dictionary as a simple argument or alternatively all keyword
        arguments will become (key, value) pairs.

        >>> meta1 = Meta({'a':1, 'b':2})
        >>> meta2 = Meta(a=1, b=2)
        """
        self.update(dict(*args, **kwargs))

    def __repr__(self):
        items = (f'{k}={v!r}' for k, v in self.__dict__.items())
        return '{}({})'.format(type(self).__name__, ', '.join(items))

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        if (isinstance(value, collections.abc.Mapping) and not
                isinstance(value, Attr)):
            self.__dict__[key] = Attr(value)
        else:
            self.__dict__[key] = value

    def __delitem__(self, name):
        del self.__dict__[name]

    def __getattr__(self, name):
        try:
            return self.__getitem__(name)
        except KeyError as e:
            raise AttributeError(e.args[0])

    __setattr__ = __setitem__
    __delattr__ = __delitem__

    def copy(self):
        return copy.deepcopy(self)

    def update(self, adict={}):
        for (key, value) in adict.items():
            self.__setitem__(key, value)

    def __iter__(self):
        return iter(self.__dict__)

    def __len__(self):
        return len(self.__dict__)


class Meta(Attr):
    """
    A class representing sequence metadata
    """
    def __str__(self):
        return self.tostr()

    def _repr_pretty_(self, p, cycle):
        if cycle:
            p.text('...')
        else:
            p.text(str(self))

    def tostr(self, w=80):
        def _key2str():
            line = f'{k:>{lenkey}}: {self[k]}'
            if len(line) > w:
                line = line[:w-3] + '...'
            out.append(line + '\n')
        out = []
        keys = set(self)
        if len(keys) == 0:
            return ''
        lenkey = max(len(k) for k in keys)
        for k in ('id',):
            if k in self:
                _key2str()
                keys.discard(k)
        for k in sorted(keys - {'features'}):
            if not k.startswith('_'):
                _key2str()
        for k in sorted(keys - {'features'}):
            if k.startswith('_'):
                _key2str()
        if 'features' in self:
            out.append(f'{"features":>{lenkey}}:\n')
            out.append(str(self.features))
            # for ft in self.features:
            #     ftstr = str(ft)
            #     if len(ftstr) > w - 25:
            #         ftstr = ftstr[:w-28] + '...'
            #     l, = ft.locs
            #     out.append('{:>13} {:<10} {}\n'.format(
            #         getattr(ft, 'type', ''), getattr(ft, 'loc', ''), ftstr))
        return ''.join(out)
