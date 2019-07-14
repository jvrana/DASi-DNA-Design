from .__version__ import __version__, __title__, __authors__, __homepage__, __repo__

from lobio.models import Region, Context
import networkx as nx
from .utils import sort_with_keys, bisect_slice_between
from more_itertools import partition, flatten
from bisect import bisect_left

class Constants(object):
    PCR_PRODUCT = "PCR_PRODUCT"
    PRIMER = "PRIMER"

    COLOR = "edge_type"
    RED = "molecule"
    BLUE = "assembly"


class AlignmentException(Exception):
    pass


class Alignment(object):
    def __init__(self, query_region, subject_region, type, query_key, subject_key):
        self.query_region = query_region
        self.subject_region = subject_region
        self.validate()
        self.type = type
        self.query_key = query_key
        self.subject_key = subject_key

    def validate(self):
        if not len(self.query_region) == len(self.subject_region):
            raise AlignmentException(
                "Regions must have the same size: {} vs {}".format(
                    len(self.query_region), len(self.subject_region)
                )
            )

    def sub_region(self, qstart, qend, type=None):
        query_copy = self.query_region.sub_region(qstart, qend)
        subject_copy = self.subject_region.copy()
        delta_s = qstart - self.query_region.start
        delta_e = self.query_region.end - qend
        subject_copy.extend_left_end(-delta_s)
        subject_copy.extend_right_end(-delta_e)
        if type is None:
            type = self.type
        self.validate()
        return self.__class__(
            query_region=query_copy,
            subject_region=subject_copy,
            type=type,
            query_key=self.query_key,
            subject_key=self.subject_key
        )

    @classmethod
    def new_pcr_product(cls):
        return cls(Constants.PCR_PRODUCT, "PCR amplify")

    def __str__(self):
        return "<{} {} {} {}>".format(
            self.__class__.__name__, self.type, self.query_region, self.subject_region
        )


class AlignmentGroup(object):

    def __init__(self, query_region, alignments):
        self.query_region = query_region
        self.alignments = alignments

    def sub_region(self, qstart, qend, type):

        return self.__class__(
            query_region=self.query_region.sub_region(qstart, qend),
            alignments=[a.sub_region(qstart, qend, type) for a in self.alignments]
        )


def blast_to_region(query_or_subject, seqdb):
    data = query_or_subject
    l = data["length"]
    record = seqdb[data["origin_key"]]
    l = len(record)
    # if data['circular']:
    #     assert l % 2 == 0
    #     l = int(l / 2.0)
    s, e = data["start"], data["end"]
    if data["strand"] == -1:
        s, e = e, s
    region = Region(
        s,
        e,
        direction=data["strand"],
        context=Context(l, data["circular"], start_index=1),
        name="{}: {}".format(record.id, record.name),
    )
    return region


class AlignmentContainer(object):

    valid_types = [Constants.PCR_PRODUCT, Constants.PRIMER]

    def __init__(self, seqdb):
        self.alignments = []
        self.seqdb = seqdb

    def load_blast_json(self, data, type):
        assert type in self.valid_types
        for d in data:
            query_region = blast_to_region(d["query"], self.seqdb)
            subject_region = blast_to_region(d["subject"], self.seqdb)
            query_key = d["query"]["origin_key"]
            subject_key = d["subject"]["origin_key"]

            alignment = Alignment(
                query_region,
                subject_region,
                type=type,
                query_key=query_key,
                subject_key=subject_key,
            )
            self.alignments.append(alignment)

    # TODO: why is there an overhang for some primers? These are supposed to be perfect alignments.
    def find_primer_pairs(self):
        pcr_products = self.groups_by_type[Constants.PCR_PRODUCT]
        primers = list(flatten([g.alignments for g in self.groups_by_type[Constants.PRIMER]]))
        rev, fwd = partition(lambda p: p.subject_region.direction == 1, primers)
        rev, rev_keys = sort_with_keys(rev, key=lambda p: p.query_region.start)
        fwd, fwd_keys = sort_with_keys(fwd, key=lambda p: p.query_region.end)

        pairs = []

        for g in pcr_products:
            fwd_bind = bisect_slice_between(fwd, fwd_keys, g.query_region.start+10, g.query_region.end)
            rev_bind = bisect_slice_between(rev, rev_keys, g.query_region.start, g.query_region.end-10)
            rkeys = [r.query_region.start for r in rev_bind]
            for f in fwd_bind:
                i = bisect_left(rkeys, f.query_region.start)
                for r in rev_bind[i:]:
                    primer_group = g.sub_region(
                        f.query_region.start,
                        r.query_region.end,
                        Constants.PCR_PRODUCT
                    )
                    pairs += primer_group.alignments
        return pairs

    def expand_pcr_products(self):
        pcr_products = self.groups_by_type[Constants.PCR_PRODUCT]
        _, product_keys = sort_with_keys(pcr_products, key=lambda x: x.query_region.start)
        for g in pcr_products:
            i = bisect_left(product_keys, g.query_region.start)
            arr, keys = sort_with_keys(pcr_products[i:], key=lambda x: x.query_region.start)
            overlapping = bisect_slice_between(arr, keys, g.query_region.start, g.query_region.end)

            for og in overlapping:
                if og is not g:
                    if og.query_region.start - g.query_region.start > 20:
                        ag1 = g.sub_region(g.query_region.start, og.query_region.start, Constants.PCR_PRODUCT)
                        self.alignments += ag1.alignments
                    if g.query_region.end - og.query_region.start > 20:
                        ag2 = g.sub_region(og.query_region.start, g.query_region.end, Constants.PCR_PRODUCT)
                        self.alignments += ag2.alignments

    def build_assembly_graph(self):
        G = nx.DiGraph()

        pcr_products = self.groups_by_type[Constants.PCR_PRODUCT]
        print("Number of alignments: {}".format(len(self.alignments)))
        pcr_products = self.groups_by_type[Constants.PCR_PRODUCT]
        for g in pcr_products:
            G.add_edge(
                g.query_region.start,
                g.query_region.end,
                **{Constants.COLOR: Constants.RED}
            )
        pairs = self.find_primer_pairs()

        print("Number of pairs: {}".format(len(pairs)))

        self.expand_pcr_products()
        print("Number of new alignments: {}".format(len(self.alignments)))
        return G

    @staticmethod
    def alignment_hash(a):
        return (
            a.query_region.start,
            a.query_region.end,
            a.query_region.direction,
            a.type,
        )

    @classmethod
    def group(cls, alignments):
        grouped = {}
        for a in alignments:
            grouped.setdefault(cls.alignment_hash(a), list()).append(a)
        alignment_groups = []
        for group in grouped.values():
            alignment_groups.append(AlignmentGroup(group[0].query_region, group))
        return alignment_groups

    @property
    def types(self):
        return tuple(self.valid_types)

    @property
    def groups_by_type(self):
        d = {}
        for t in self.valid_types:
            d[t] = []
        for a in self.alignments:
            d[a.type].append(a)
        return {t: self.group(v) for t, v in d.items()}

    @property
    def alignment_groups(self):
        return self.group(self.alignments)
