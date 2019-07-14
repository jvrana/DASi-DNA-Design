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
        if not len(query_region) == len(subject_region):
            raise AlignmentException(
                "Regions must have the same size: {} vs {}".format(
                    len(query_region), len(subject_region)
                )
            )
        self.type = type
        self.query_key = query_key
        self.subject_key = subject_key

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

    def build_assembly_graph(self):
        G = nx.DiGraph()

        pcr_products = self.groups_by_type[Constants.PCR_PRODUCT]
        primers = list(flatten([g.alignments for g in self.groups_by_type[Constants.PRIMER]]))
        for g in pcr_products:
            G.add_edge(
                g.query_region.start,
                g.query_region.end,
                **{Constants.COLOR: Constants.RED}
            )

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
                for r in rev_bind[:]:
                    pairs.append((f, r))

        print(len(pairs))
        return G


    # def index_hash_table(self, alignments):
    #     by_start = sorted(pcr_products, key=lambda a: a.query_region.start)
    #     by_end = sorted(pcr_product, key=lambda a: a.query_region.end)

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
