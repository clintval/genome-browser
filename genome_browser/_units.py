import csv
import gzip
import io
import re

__all__ = [
    'Graph',
    'Interval',
    'Exon',
    'Gene',
    'RefGene']


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
        assert end - start > 0, 'Exclusive end must be greater than start'
        assert strand in ('.', '+', '-'), 'Strand must be ".", "+", "-" only'

        self.chrom = str(chrom)
        self.start = start
        self.end = end
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


class Exon(Interval):
    def __init__(
        self,
        chrom,
        start,
        end,
        strand='.',
        rank=None,
        frame_offset=-1,
        name=None,
        **metadata
    ):
        super().__init__(
            chrom,
            start,
            end,
            strand=strand,
            name=name,
            metadata=metadata)

        if (
            frame_offset is not None and
            not (isinstance(frame_offset, int) and
                 frame_offset in range(-1, 3))
        ):
            raise ValueError(
                'Frame offset must be None or integer and in set [-1, 2]')

        if not isinstance(rank, int) and rank > 0:
            raise ValueError('Rank must be a positive integer!')

        self.rank = rank
        self.frame_offset = None if frame_offset == -1 else frame_offset

    @property
    def in_frame(self):
        return self.frame_offset == 0

    @property
    def has_frame(self):
        return self.frame_offset is not None

    def __repr__(self):
        return (
            'Exon('
            f'"{self.chrom}", '
            f'{self.start}, '
            f'{self.end}, '
            f'"{self.strand}", '
            f'rank={self.rank}, '
            f'frame_offset={self.frame_offset})')


class Gene(Interval):
    def __init__(
        self,
        chrom,
        start,
        end,
        strand='.',
        name=None,
        id=None,
        coding_start=None,
        coding_end=None,
        score=None,
        coding_start_status=None,
        coding_end_status=None,
        **metadata
    ):
        super().__init__(
            chrom,
            start,
            end,
            strand=strand,
            name=name,
            metadata=metadata)

        self.id = id
        self.transcript_start = start
        self.transcript_stop = end
        self.coding_start = coding_start
        self.coding_end = coding_end
        self.coding_start_status = coding_start_status
        self.coding_end_status = coding_end_status

        self._exons = []

    @property
    def num_exons(self):
        return len(self.exons)

    @property
    def exons(self):
        return self._exons

    @exons.getter
    def exons(self):
        return sorted(self._exons)

    def __repr__(self):
        return (
            'Gene('
            f'"{self.chrom}", '
            f'{self.start}, '
            f'{self.end}, '
            f'"{self.strand}", '
            f'name="{self.name}", '
            f'id="{self.id}")')


class RefGene(object):
    def __init__(self, gzip_file):
        self.gzip_file = gzip_file

    @staticmethod
    def _line_to_gene(line):
        (_, name, chrom, strand, transcript_start, transcript_stop,
         coding_start, coding_end, num_exons, exon_starts, exon_ends,
         score, alt_name, coding_start_status, coding_end_status,
         exon_frames) = line

        gene = Gene(
            chrom,
            int(transcript_start),
            int(transcript_stop),
            strand,
            name=alt_name,
            id=name,
            coding_start=int(coding_start),
            coding_end=int(coding_end),
            score=score,
            coding_start_status=coding_start_status,
            coding_end_status=coding_end_status)

        exon_starts = exon_starts.split(',')
        exon_ends = exon_ends.split(',')
        exon_frames = exon_frames.split(',')

        if strand == '-':
            exon_ranks = range(int(num_exons), 0, -1)
        else:
            exon_ranks = range(1, int(num_exons) + 1)

        for start, end, frame_offset, rank in zip(
            exon_starts, exon_ends, exon_frames, exon_ranks
        ):
            if any(_ == '' for _ in (start, end, frame_offset)):
                continue

            exon = Exon(
                chrom=chrom,
                start=int(start),
                end=int(end),
                strand=strand,
                rank=rank,
                frame_offset=int(frame_offset))

            gene._exons.append(exon)

        return gene

    def gene_by_id(self, id):
        for gene in self:
            if gene.id == id:
                return gene

    def genes_by_id_pattern(self, id, flags=re.IGNORECASE):
        pattern = re.compile(id, flags)

        for gene in self:
            if pattern.match(gene.id):
                yield gene

    def gene_by_name(self, name):
        for gene in self:
            if gene.name == name:
                return gene

    def genes_by_name_pattern(self, name, flags=re.IGNORECASE):
        pattern = re.compile(name, flags)

        for gene in self:
            if pattern.match(gene.name):
                yield gene

    def __iter__(self):
        handle = gzip.open(self.gzip_file, 'rb')
        self._reader = csv.reader(
            io.TextIOWrapper(handle),
            delimiter='\t')
        return self

    def __next__(self):
        return self._line_to_gene(next(self._reader))

    def __repr__(self):
        return f'RefSeq("{self.gzip_file}")'
