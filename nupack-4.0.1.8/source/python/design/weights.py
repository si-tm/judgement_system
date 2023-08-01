import pandas, numpy
from ..utility import printable, long_output
from ..core import TargetComplex, TargetTube

################################################################################

class Weight:
    '''Individual weighting factor'''
    def __init__(self, tube, complex, strand, domain, weight):
        pass

################################################################################

class ReversedComplex:
    '''Helper class for complex weights'''

    def strands(self):
        pass

    def domains(self):
        pass

################################################################################

def _default_format():
    return '{}'

def _get_name(x):
    return x.name

def style(table, names, html, **formats):
    '''Format each column by the given formats'''
    if names is not None:
        table = table[table.columns[:len(names)]]
        formats.update({k: formats[v] for k, v in zip(names, table.columns) if v in formats})
        table.columns = names

    if html:
        return table.to_html(formatters=formats, index=False, na_rep='')
    else:
        return table.to_string(formatters=formats, index=False, na_rep='')

################################################################################

class Array:
    '''
    Class representing essemtially a pandas DataFrame with a MultiIndex key
    of [tube, complex, strand, domain]. Access that DataFrame as `.named()`.

    The actual table is not held as a MultiIndex for a simpler implementation.
    '''
    table: pandas.DataFrame


    def __init__(self, tubes_or_complexes):
        '''Initialize from a list of tubes'''
        if all(isinstance(t, TargetTube) for t in tubes_or_complexes):
            rows = sorted(set((d, s, c, t) for t in tubes_or_complexes
                for c in t.on_targets for s in c.strands for d in s.domains))
            self.axes = ('domain', 'strand', 'complex', 'tube')
        elif all(isinstance(t, TargetComplex) for t in tubes_or_complexes):
            rows = sorted(set((d, s, c) for c in tubes_or_complexes for s in c.strands for d in s.domains))
            self.axes = ('domain', 'strand', 'complex')
        else:
            raise TypeError('Expected list of TargetTube or list of TargetComplex')
        self.table = pandas.DataFrame(rows, columns=self.axes)

    def replace(self, rules):
        self.table = self.table.applymap(lambda x: (
            x.substitute(rules) if hasattr(x, 'substitute') else x))

    def named(self, upper=False):
        '''Return an equivalent DataFrame using a MultiIndex'''
        df = self.table.copy()
        for a in self.axes:
            df[a] = df[a].map(_get_name)
        if upper:
            df.columns = [s[0].upper() + s[1:] for s in df.columns]
        return df

    def __str__(self):
        with long_output():
            return style(self.named(True), None, False, weight='{:#.3g}'.format)

    def _repr_html_(self):
        with long_output():
            return style(self.named(True), None, True, weight='{:#.3g}'.format)

    def __getattr__(self, key):
        if key.startswith('_'):
            raise AttributeError(key)
        return getattr(self.table, key)

    def mask(self, key):
        '''Return rows which match the specified key'''
        if isinstance(key, tuple):
            key = key + (len(self.axes) - len(key)) * (slice(None),)
        else:
            key = (key,) + (len(self.axes) - 1) * (slice(None),)
        mask = numpy.full(len(self.table), True)
        for a, k in zip(self.axes, key):
            if isinstance(k, slice):
                bools = numpy.full(len(self.table), False)
                bools[k] = True
                mask &= bools
            elif isinstance(k, str):
                mask &= [f.name == k for f in self.table[a]]
            else:
                mask &= (self.table[a] == k)
        if not numpy.any(mask):
            raise KeyError(key)
        return mask

################################################################################

@printable
class Weights(Array):
    '''
    Class representing objective weights on domains, strands, complexes, and tubes
    '''

    def __init__(self, tubes_or_complexes):
        '''Initialize to 1s from a list of tubes or list of complexes'''
        super().__init__(tubes_or_complexes)
        self.table['weight'] = 1.0

    def __getitem__(self, key):
        '''Get weights matching a specified key'''
        return self.table.loc[self.mask(key), 'weight']

    def __setitem__(self, key, value):
        '''Set weights matching a key to value'''
    #     assert numpy.all(value >= 0), 'Weights must be non-negative'
        self.table.loc[self.mask(key), 'weight'] = value

    def factors(self):
        return [Weight(r['tube'].name, r['complex'].name, r['strand'].name, r['domain'].name, r['weight']) 
                for _, r in self.table[self.table.weight != 1].iterrows()]
