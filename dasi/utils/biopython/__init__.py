"""biopython.py.

Helper functions for BioPython
"""
import hashlib
import itertools
import random
import re
import tempfile
from copy import deepcopy
from itertools import chain
from typing import List
from typing import Tuple
from typing import Union
from uuid import uuid4

import networkx as nx
from Bio import Restriction
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import ExactPosition
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from primer3plus.design import primer3

from .record_to_json import seqrecord_to_json
from dasi.utils.biopython.sequence import anneal
from dasi.utils.biopython.sequence import random_sequence
from dasi.utils.biopython.sequence import rc
from dasi.utils.region import Span

FWD_COLOR = "ApEinfo_fwdcolor"
REV_COLOR = "ApEinfo_revcolor"
CIRCULAR = "circular"
LINEAR = "lienar"
TOPOLOGY = "topology"
DEFAULT_ANNOTATION_TYPE = "misc_feature"


def sort_cycle(arr, key=None):
    """Sort a cyclic array, maintaining order."""
    if key is None:
        arr_with_i = sorted([(x, i) for i, x in enumerate(arr)])
    else:
        arr_with_i = sorted([(key(x), i) for i, x in enumerate(arr)])
    i = arr_with_i[0][1]
    return arr[i:] + arr[:i]


def random_color():
    """Make a random color."""
    random_number = random.randint(0, 16777215)
    hex_number = str(hex(random_number))[2:]
    if len(hex_number) < 6:
        hex_number = "0" * (6 - len(hex_number)) + hex_number
    return "#" + hex_number


def random_slices(mn: int, mx: int, total_length: int):
    """Yield random slices whose lengths sum to the provided total length.

    :param mn: minimum slice size
    :param mx: maximum slice size
    :param total_length: total length of the slices
    :return:
    """
    n = total_length
    j = 0
    while j < n:
        i = j
        j += random.randint(mn, mx)
        if j >= n:
            j = n
        yield (i, j)


def format_float(a, places=2):
    return "%." + str(places) + "f" % round(a, places)


def get_feature_name(feature: SeqFeature):
    return get_feature_qualifiers(feature, "label")


def set_feature_name(feature: SeqFeature, name: str):
    set_feature_qualifier(feature, "label", name)
    return feature


def set_feature_qualifier(feature: SeqFeature, key: str, value):
    feature.qualifiers[key] = [value]


def get_feature_qualifiers(feature: SeqFeature, key: str):
    return feature.qualifiers.get(key, list())[0]


def set_feature_color(feature: SeqFeature, color: str = None, rev_color: str = None):
    if color is None:
        color = random_color()

    if not rev_color:
        rev_color = random_color()

    set_feature_qualifier(feature, FWD_COLOR, color)
    set_feature_qualifier(feature, REV_COLOR, rev_color)


def new_location(i: int, j: int, strand: int) -> FeatureLocation:
    return FeatureLocation(ExactPosition(i), ExactPosition(j), strand=strand)


def new_compound_location(
    indices: List[Union[Tuple[int, int], Tuple[int, int, int]]], strand: int
) -> CompoundLocation:
    locations = []
    for index in indices:
        if not isinstance(index, Tuple):
            raise ValueError(
                "Expects a tuple of integers size 2 or 3, not a {}".format(
                    indices.__class__
                )
            )
        if not len(index) in [2, 3]:
            raise ValueError("Expects a tuple of integers of size 2 or 3")
        if len(index) == 2:
            i, j = index
            s = strand
        elif len(index) == 3:
            i, j, s = index
        else:
            raise ValueError("Must be tuple of 2 or 3 integers")
        if not isinstance(i, int) or not isinstance(j, int) or not isinstance(s, int):
            raise ValueError(
                "Expects a tuple of integers of size 2 or 3. Found {}".format(index)
            )
        locations.append(FeatureLocation(ExactPosition(i), ExactPosition(j), strand=s))
    return CompoundLocation(locations)


def clean_features(record: SeqRecord):
    """Remove redundant features."""
    features_by_hash = {}
    for f in record.features:
        fhash = "*".join([str(f.qualifiers["label"]), f.type, str(f.location)])
        features_by_hash[fhash] = f
    record.features = list(features_by_hash.values())
    return record


def new_compound_feature(
    name: str,
    indices: List[Union[Tuple[int, int], Tuple[int, int, int]]],
    strand: int,
    feature_type: str = None,
):
    location = new_compound_location(indices, strand=strand)
    return new_feature(name, location=location, feature_type=feature_type)


def new_feature(
    name: str,
    location: Union[None, CompoundLocation, FeatureLocation],
    color: str = None,
    feature_type: str = None,
):
    feature = SeqFeature(location=location, qualifiers={}, type=feature_type)
    set_feature_name(feature, name)
    set_feature_color(feature, color=color)
    return feature


def _annotate_feature(
    length: int,
    name: str,
    i: int = None,
    j: int = None,
    cyclic: bool = False,
    feature_type: str = None,
):
    if i is None:
        i = 0
    if j is None:
        j = length

    if cyclic and (j > length or (i > j)):
        if j > length:
            j = j - length
        feature = new_compound_feature(
            name=name,
            indices=[(i, length), (0, j)],
            strand=1,
            feature_type=feature_type,
        )
    else:
        feature = new_feature(
            name=name,
            location=FeatureLocation(ExactPosition(i), ExactPosition(j), strand=1),
            feature_type=feature_type,
        )
    return feature


def annotate(
    record: SeqRecord,
    name: str,
    i: int = None,
    j: int = None,
    cyclic: bool = False,
    annotation_type: str = None,
):
    """Annotate a SeqRecord."""
    if not name:
        raise ValueError("Cannot annotate record with no name.")
    if annotation_type is None:
        annotation_type = DEFAULT_ANNOTATION_TYPE
    feature = _annotate_feature(
        len(record.seq), name, i, j, cyclic, feature_type=annotation_type
    )
    record.features.append(feature)
    return feature


def set_topology(records: List[SeqRecord], topology: str):
    assert topology in [CIRCULAR, LINEAR]
    for r in records:
        r.annotations[TOPOLOGY] = topology


def make_cyclic(records: List[SeqRecord]):
    if isinstance(records, SeqRecord):
        records = [records]
    for r in records:
        r.annotations[TOPOLOGY] = CIRCULAR


def make_linear(records: List[SeqRecord]):
    if isinstance(records, SeqRecord):
        records = [records]
    for r in records:
        r.annotations[TOPOLOGY] = LINEAR


def is_circular(record: SeqRecord):
    return record.annotations.get(TOPOLOGY, False) == CIRCULAR


def is_linear(record: SeqRecord):
    return not is_circular(record)


def new_sequence(
    seqstr: str,
    name: str = None,
    auto_annotate: bool = False,
    cyclic: bool = False,
    annotation_type: str = None,
):
    record = SeqRecord(Seq(seqstr))
    if auto_annotate:
        annotate(record, name, annotation_type=annotation_type)
    if name:
        record.name = name
        record.id = name
        record.description = name

    if cyclic:
        make_cyclic([record])
    else:
        make_linear([record])
    return record


# def slice_feature(feature: SeqFeature, i: int = None, j: int = None):
#     if isinstance(feature)


def _format_slice(slc: slice):
    return "[{}:{}]".format(slc.start, slc.stop)


def slice_with_features(seq: SeqRecord, slc: slice):
    i = slc.start
    j = slc.stop
    if slc.step:
        raise ValueError("Step is not supported.")
    # new_features = []
    # for feature in seq.features:
    #     name = get_feature_name(feature)
    #     new_parts = []
    #     for part in feature.location.parts:
    #         if i in part and j in part:
    #             new_parts.append(new_location(i, j, part.strand))
    #         elif i in part:
    #             new_parts.append(new_location(i, part.end.position, part.strand))
    #         elif j in part:
    #             new_parts.append(new_location(part.start.position, j, part.strand))
    #     if len(new_parts) == 1:
    #         new_feature = deepcopy(feature)
    #         new_feature.location = new_parts[0]
    #         set_feature_name(new_feature, name + _format_slice(slc))
    #         new_features.append(new_feature)
    #     elif len(new_parts) > 2:
    #         new_feature = deepcopy(feature)
    #         new_feature.location = CompoundLocation(new_parts)
    #         set_feature_name(new_feature, name + _format_slice(slc))
    #         new_features.append(new_feature)
    new_seq = seq[i:j]
    # new_seq.features += new_features
    return new_seq


def random_record(
    length: int, name: str = None, auto_annotate: bool = True, cyclic: bool = False
):
    if name is None:
        name = str(uuid4())
    return new_sequence(
        random_sequence(length), name=name, auto_annotate=auto_annotate, cyclic=cyclic
    )


def randomly_annotate(
    rec: SeqRecord,
    feature_length_range: Tuple[int, int],
    feature_name_list: List[str] = None,
) -> SeqRecord:
    for i, j in random_slices(
        feature_length_range[0], feature_length_range[1], len(rec.seq)
    ):
        if not feature_name_list:
            random_name = str(uuid4())[-4:]
        else:
            random_name = random.sample(feature_name_list, k=1)[0]
        annotate(rec, random_name, i, j)
    return rec


def pcr_amplify(
    fwd: SeqRecord,
    rev: SeqRecord,
    template: SeqRecord,
    cyclic: bool,
    cyclic_buffer: int = 100,
    name: str = None,
    return_matches: bool = False,
):
    original_template = template
    if cyclic:
        template = template + template + template[:cyclic_buffer]
    fwd_matches, rev_matches = anneal(str(template.seq), [str(fwd.seq), str(rev.seq)])

    products_by_sequence = {}
    for f, r in itertools.product(fwd_matches, rev_matches):
        i = f["top_strand_slice"][0]
        j = r["top_strand_slice"][1]

        try:
            span = Span(
                i,
                j,
                length=len(original_template.seq),
                cyclic=cyclic,
                ignore_wrap=True,
            )
        except IndexError as e:
            if not cyclic:
                continue
            else:
                raise e

        f_rec = new_sequence(
            f["overhang"], name=fwd.name + "_overhang", auto_annotate=True
        )
        r_rec = new_sequence(
            rc(r["overhang"]), name=rev.name + "_overhang", auto_annotate=True
        )

        product = span.get_slice(original_template)
        # annotate(template_chunk, name="source: {}".format(template_name[:40]))
        if len(f_rec.seq):
            product = f_rec + product
        if len(r_rec.seq):
            product = product + r_rec

        if not name:
            name = original_template.name or original_template.id
            name += "[{}:{}]".format(span.a, span.b)
        # annotate(product, name=name)
        if len(product) <= len(original_template.seq):
            products_by_sequence[str(product.seq)] = (product, span.a, span.b)
    product_list = list(products_by_sequence.values())
    product_list.sort(key=lambda x: (x[1], len(product)))
    if return_matches:
        return product_list, fwd_matches, rev_matches
    return product_list


def load_glob(*paths: Tuple[str, ...], format: str = None):
    path_iter = chain(*paths)
    records = []
    for path in path_iter:
        records += SeqIO.parse(path, format=format)
    return records


def load_fasta_glob(path: str) -> List[SeqRecord]:
    return load_glob(path, format="fasta")


def load_genbank_glob(path: str) -> List[SeqRecord]:
    return load_glob(path, format="genbank")


def write_tmp_records(records: List[SeqRecord], format: str) -> List[SeqRecord]:
    fd, tmp_path_handle = tempfile.mkstemp(suffix="." + format)
    SeqIO.write(records, tmp_path_handle, format=format)
    return tmp_path_handle


class GibsonAssembler:
    @staticmethod
    def make_hash(s: str):
        return hashlib.sha1(s.encode("utf-8")).hexdigest()

    @classmethod
    def rnode(cls, record: SeqRecord):
        return cls.make_hash(str(record.name + record.seq))

    @classmethod
    def add_edge(
        cls, g: nx.DiGraph, r1: SeqRecord, r2: SeqRecord, tm: float, match: dict
    ):
        n1 = cls.rnode(r1)
        n2 = cls.rnode(r2)
        if n1 not in g:
            g.add_node(n1, record=r1)
        if n2 not in g:
            g.add_node(n2, record=r2)
        g.add_edge(n1, n2, tm=tm, match=match)

    @classmethod
    def interaction_graph(cls, records: List[SeqRecord]) -> nx.DiGraph:
        g = nx.DiGraph()
        pairs = []
        for r1, r2 in itertools.product(records, records):
            if r1 is not r2:
                pairs.append((r1, r2))
                if r2 is not records[0]:
                    pairs.append((r1, r2.reverse_complement(name=r2.name + "_rc")))

        for r1, r2 in pairs:
            fwd, rev = anneal(str(r1.seq), [str(r2.seq)])
            for f in fwd:
                if f["top_strand_slice"][0] == 0:
                    anneal_seq = f["anneal"]
                    tm = primer3.calcTm(anneal_seq[-60:])
                    cls.add_edge(g, r2, r1, tm=tm, match=f)
        return g

    @classmethod
    def _collect_features(cls, records: List[SeqRecord]):
        seq_to_features = {}
        for record in records:
            for feature in record.features:
                s = str(feature.extract(record.seq)).upper()
                if s not in seq_to_features:
                    seq_to_features[s] = feature
        return seq_to_features

    @classmethod
    def _restore_features(cls, seq_to_features, record, cyclic: bool):
        for seq in seq_to_features:
            record_seq = str(record.seq)
            if cyclic:
                record_seq += record_seq
            matches = list(re.finditer(seq, record_seq, re.IGNORECASE))
            if len(matches) >= 1:
                feature = seq_to_features[seq]
                new_feature = deepcopy(feature)
                for match in matches:
                    span = match.span()

                    if span[0] < len(record.seq):
                        if span[1] >= len(record.seq):
                            f = _annotate_feature(
                                len(record.seq), "", span[0], span[1], cyclic=True
                            )
                            new_feature.location = f.location
                        else:
                            new_feature.location = new_location(
                                span[0], span[1], feature.location.strand
                            )
                    record.features.append(new_feature)

    @classmethod
    def assemble_record_from_cycle(
        cls,
        cycle: List[str],
        graph: nx.DiGraph,
        annotate_sources: bool = False,
        annotate_junctions: bool = False,
    ) -> SeqRecord:

        source_indices = []
        i = 0

        records = [graph.nodes[n]["record"] for n in cycle]
        stored_features = cls._collect_features(records)

        record = SeqRecord(Seq(""))

        c1 = cycle[-1:] + cycle[:-1]
        c2 = cycle
        c3 = cycle[1:] + cycle[:1]

        for n1, n2, n3 in zip(c1, c2, c3):
            r1, r2, r3 = [graph.nodes[n]["record"] for n in [n1, n2, n3]]

            edata12 = graph.get_edge_data(n1, n2)
            edata23 = graph.get_edge_data(n2, n3)

            s12 = edata12["match"]["top_strand_slice"]
            s23 = edata23["match"]["top_strand_slice"]

            o12 = r2[s12[0] : s12[1]]
            o23 = r3[s23[0] : s23[1]]

            trunc_r2 = slice_with_features(r2, slice(len(o12), -len(o23)))

            if annotate_junctions:
                annotate(o12, "junction: tm {}".format(format_float(edata12["tm"])))
            record += o12
            record += trunc_r2

            source_indices.append((i, i + len(r2)))
            i = i + len(r2) - len(o23)

        for i, (a, b) in enumerate(source_indices):
            if annotate_sources:
                annotate(
                    record,
                    "source: fragment {} ({})".format(i, r2.name),
                    a,
                    b,
                    cyclic=True,
                )

        cls._restore_features(stored_features, record, cyclic=True)
        clean_features(record)
        return record

    @classmethod
    def make_cyclic_assemblies(
        cls,
        records: List[SeqRecord],
        annotate_sources: bool = False,
        annotate_junctions: bool = False,
    ) -> List[SeqRecord]:
        g = cls.interaction_graph(records)

        node_rank = {}
        for i, r in enumerate(records):
            node_rank[cls.rnode(r)] = i

        cycles = list(nx.simple_cycles(g))
        cycles = [sort_cycle(cycle, key=lambda n: node_rank[n]) for cycle in cycles]

        assembled_records = []
        for cycle in cycles:
            record = cls.assemble_record_from_cycle(
                cycle,
                g,
                annotate_sources=annotate_sources,
                annotate_junctions=annotate_junctions,
            )
            assembled_records.append(record)
        make_cyclic(assembled_records)
        return assembled_records


def digest(record: SeqRecord, restriction_enzyme: str):
    """Digest a SeqRecord with an enzyme.

    :param record: record to digest
    :param restriction_enzyme: name of restriction enzyme
    :return: list of records
    """
    enzyme = getattr(Restriction, restriction_enzyme)
    linear = is_linear(record)
    indices = enzyme.search(record.seq, linear=linear)
    indices = [i - 1 for i in indices]
    pairs = zip(indices[:-1], indices[1:])
    records = []
    for i1, i2 in pairs:
        records.append(record[i1:i2])
    if not linear:
        records.append(record[indices[-1] :] + record[: indices[0]])
    return records


def random_record_from_library(
    records: List[SeqRecord],
    circular: bool,
    size_interval: Tuple[int, int] = (5000, 10000),
    max_chunks: int = None,
    chunk_size_interval: Tuple[int, int] = (100, 3000),
    random_chunk_prob_int: Tuple[float, float] = (0, 0.5),
    random_chunk_size_int: Tuple[int, int] = (100, 1000),
):
    """

    :param records:
    :param size_interval:
    :param chunk_size_interval:
    :param random_chunk_prob_interval: picks a random probability at whic
    :param circular:
    :return:
    """
    length = random.randint(*size_interval)
    new_record = new_sequence("")

    x = random.uniform(*random_chunk_prob_int)

    n_chunks = 0
    while len(new_record.seq) < length and (
        max_chunks is None or n_chunks < max_chunks
    ):
        if random.random() < x:
            rand_chunk_len = random.randint(*random_chunk_size_int)
            rand_rec = random_record(rand_chunk_len)
        else:
            rand_chunk_len = random.randint(*chunk_size_interval)
            rand_rec = random.choice(records)
            if rand_chunk_len > len(rand_rec.seq):
                continue
                # rand_chunk_len = len(rand_rec.seq)

        i = random.randint(0, len(rand_rec.seq) - rand_chunk_len)
        rand_chunk = rand_rec[i : i + rand_chunk_len]
        if random.random() < 0.5:
            rand_chunk = rand_chunk.reverse_complement()

        new_record += rand_chunk
        n_chunks += 1
    new_record = new_record[:length]
    if circular:
        make_cyclic([new_record])
    else:
        make_linear([new_record])
    return new_record


make_cyclic_assemblies = GibsonAssembler.make_cyclic_assemblies
