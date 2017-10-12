import csv
import json
import requests

import matplotlib.pyplot as plt

from matplotlib.patches import Rectangle, Shadow
from palettable.cartocolors.qualitative import Prism_10
# from palettable.colorbrewer.qualitative import Paired_12

from genome_browser._util import ax_off, despine, string_as_temporary_file

__all__ = [
    'get_uniprot_id_from_gene_symbol',
    'get_pfam_entry_from_uniprot_id',
    'get_pfam_entry_graphic',
    'get_features_from_uniprot_id',
    'plot_gene_features_from_seq_record']


FACECOLOR = '0.978'
LABEL_SPACE = 5
TEXT_ZORDER = 20
SHADOW_PROPS = {'alpha': 0.075, 'color': 'k', 'ec': 'none'}

COLOR_MAP = {
    'dna binding': Prism_10.hex_colors[-2],
    'helix': Prism_10.hex_colors[1],
    'chain': '0.675',
}

ALPHA_MAP = {
    'dna binding': 1,
    'helix': 1,
    'chain': 1,
}

ZORDER_MAP = {
    'dna binding': 7,
    'helix': 6,
    'chain': 3,
}

HEIGHT_PADDING_MAP = {
    'dna binding': 0.2,
    'helix': 0,
    'chain': -0.75,
}

MOTIF_DEFINITION = {
    'disorder': 'Disordered region (Pfam/IUPred)',
    'low_complexity': 'Low complexity region (Pfam/SEG)',
    'sig_p': 'Signal peptide region (Pfam/Phobius)',
    'coiled_coil': 'Coiled-coil motif (Pfam/ncoils)',
    'transmembrane': 'Transmembrane region (Pfam/Phobius)'}


class PfamGraphicFeature(object):
    def __init__(self, feature):
        self.color      = feature.get('colour', 'grey')
        self.display    = feature.get('display', None)
        self.end        = feature.get('end', None)
        self.endstyle   = feature.get('endStyle', None)
        self.link       = feature.get('href', None)
        self.start      = feature.get('start', None)
        self.startstyle = feature.get('startStyle', None)
        self.text       = feature.get('text', '')
        self.type       = feature.get('type', None)
        self.metadata   = feature.get('metadata', {})

    def __repr__(self):
        return (
            f'PfamGraphicFeature('
            f'start={self.start} '
            f'end={self.end} '
            f'color="{self.color}" '
            f'link="{self.link}")')


class PfamGraphicResponse(object):
    def __init__(self, content):
        self.regions = []

        self.length   = content.get('length', None)
        self.markups  = content.get('markups', ())
        self.metadata = content.get('metadata', {})
        self.motifs   = content.get('motifs', ())

        for region in content.get('regions', ()):
            self.regions.append(PfamGraphicFeature(region))

    def __repr__(self):
        return (
            f'{self.__class__.__name__}(\n'
            f'  accession:    "{self.metadata.get("accession", "N/A")}"\n'
            f'  identifier:   "{self.metadata.get("identifier", "N/A")}"\n'
            f'  organism:     "{self.metadata.get("organism", "N/A")}")\n'
            f'  description:  "{self.metadata.get("description", "N/A")}"\n'
            f'  number of motifs:  {len(self.motifs)}\n'
            f'  number of regions: {len(self.regions)})')


def get_uniprot_id_from_gene_symbol(symbol, email, species='all'):
    from mygene import MyGeneInfo

    response = MyGeneInfo().query(
        symbol,
        fields='uniprot',
        species=species,
        email=email)

    for hit in response.get('hits'):
        uniprot_id = hit.get('uniprot', {}).get('Swiss-Prot', None)
        if uniprot_id is not None:
            break

    return uniprot_id


def get_pfam_entry_graphic(pfam_entry):
    with requests.Session() as session:
        response = session.get(
            f'http://pfam.xfam.org/protein/{pfam_entry}/graphic')
        content, *_ = json.loads(response.content.decode('utf-8'))

    return PfamGraphicResponse(content)


def get_pfam_entry_from_uniprot_id(uniprot_id):
    url = (
        f'http://www.uniprot.org/uniprot/?query={uniprot_id}'
        '+AND+reviewed:yes+AND+AND+database:pfam'
        '&sort=score&columns=entry+name,reviewed,genes,organism&format=tab')

    with requests.Session() as session:
        reader = csv.reader(
            session.get(url).content.decode('utf-8').splitlines(),
            delimiter='\t')

        next(reader)  # noqa
        entry_name, reviewed, genes, organism = next(reader)

    return entry_name


def get_features_from_uniprot_id(uniprot_id):
    url = f'http://www.uniprot.org/uniprot/{uniprot_id}.gff'

    with requests.Session() as session:
        content = session.get(url).content.decode('utf-8')

    def remove_tabs_from_line_ends(content):
        formatted = []
        for line in content.split('\n'):
            line = line.strip()
            formatted.append(line)
        return '\n'.join(formatted)

    try:
        from BCBio import GFF
        handle = string_as_temporary_file(remove_tabs_from_line_ends(content))
        records = list(GFF.parse(handle.name))
        return records
    except ValueError:  # Gosh I hate GFF parsers failing!
        return content
    except ImportError:
        return content


def plot_gene_features_from_seq_record(seq_record, ax=None):
    ax = ax or plt.gca()

    xticks = []
    x_min, x_max = float('inf'), float('-inf')

    for feature in sorted(
        seq_record.features,
        key=lambda _: _.location.end.position,
        reverse=True
    ):
        start = feature.location.start.position
        end = feature.location.end.position

        x_min, x_max = min(start, x_min), max(end, x_max)

        color = COLOR_MAP.get(feature.type.lower(), '0.7')
        zorder = ZORDER_MAP.get(feature.type.lower(), 1)
        alpha = ALPHA_MAP.get(feature.type.lower(), 1)
        height_padding = HEIGHT_PADDING_MAP.get(feature.type.lower(), 0)

        rectangle = Rectangle(
            xy=(start, -0.5 - (height_padding / 2)),
            width=end - start,
            height=1 + height_padding,
            color=color,
            zorder=zorder,
            alpha=alpha)
        ax.add_patch(rectangle)

        shadow = Shadow(
            rectangle,
            ox=0.06,
            oy=-0.007,
            props=SHADOW_PROPS)
        shadow.set_zorder(1)

        # Add shadows to all features except chains.
        if feature.type.lower() != 'chain':
            ax.add_patch(shadow)

        # Add tick marks. Starts are guaranteed, ends only appear if they are
        # not within LABEL_SPACE of another start.
        if feature.type.lower() in ('helix', 'chain'):
            if len(xticks) != 0 and xticks[-1] - end < LABEL_SPACE:
                xticks.append(start)
            else:
                xticks.extend((end, start))

        if feature.type.lower() == 'dna binding':
            ax.annotate(
                s=feature.type,
                xy=(start + (end - start) / 2, 0),
                color='w',
                ha='center',
                va='center',
                fontsize=14,
                zorder=TEXT_ZORDER)

    ax.set_xticks(sorted(set(xticks)))
    ax.set_xticklabels(sorted(set(xticks)), fontsize=9)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(-0.75, 0.75)

    despine(ax_off(ax, 'y'), bottom=False)

    return ax
