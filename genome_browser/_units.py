__all__ = [
    'Graph',
    'Interval']


class Graph(object):
    def __init__(
        self,
        x,
        y,
        fmt='interpolate',
        fill=False,
        color='0.2',
        alpha=1
    ):
        self.x = x
        self.y = y
        self.fmt = fmt
        self.fill = fill
        self.color = color
        self.alpha = alpha

    @classmethod
    def from_bedgraph(self, infile):
        ...


class Interval(object):
    def __init__(
        self,
        chrom,
        start,
        end,
        strand='.',
        name=None,
        **metadata
    ):
        assert isinstance(start, int), 'Loci must be integers'
        assert isinstance(end, int), 'Loci must be integers'
        assert strand in ('.', '+', '-'), 'Strand must be ".", "+", "-" only'

        self.chrom = str(chrom)
        self.start = min(start, end)
        self.end = max(start, end)
        self.strand = strand
        self.name = name
        self.metadata = metadata

    @property
    def sam_interval(self):
        return f'{self.chrom}:{self.start}-{self.end}'

    def get(self, item, default=None):
        temp = self.__dict__.copy()
        temp.update(self.metadata)
        return temp.get(item, default)

    def __getitem__(self, item):
        temp = self.__dict__.copy()
        temp.update(self.metadata)
        return temp.get(item, None)

    def __len__(self):
        return self.end - self.start

    def __eq__(self, other):
        return all((
            self.chrom == other.chrom,
            self.start == other.start,
            self.end == other.end,
            self.strand == other.strand))

    def __lt__(self, other):
        return all((
            self.chrom == other.chrom,
            self.start < other.start))

    def __le__(self, other):
        return all((
            self.chrom == other.chrom,
            self.start <= other.start))

    def __gt__(self, other):
        return all((
            self.chrom == other.chrom,
            self.end > other.end))

    def __ge__(self, other):
        return all((
            self.chrom == other.chrom,
            self.end >= other.end))

    def __repr__(self):
        return (
            'Interval('
            f'"{self.chrom}", '
            f'{self.start}, '
            f'{self.end}, '
            f'"{self.strand}")')

    def __str__(self):
        return self.__repr__()
