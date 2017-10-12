import csv
import json
import requests

from genome_browser._util import string_as_temporary_file

__all__ = [
    'get_uniprot_id_from_gene_symbol',
    'get_pfam_entry_from_uniprot_id',
    'get_pfam_entry_graphic',
    'get_features_from_uniprot_id']


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
