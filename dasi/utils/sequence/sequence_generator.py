"""sequence_generator."""
import random
from typing import Dict
from typing import List
from typing import Tuple
from typing import Union
from uuid import uuid4

from Bio.SeqRecord import SeqRecord
from pyblast.utils import make_circular
from pyblast.utils import make_linear

from dasi.utils import biopython


def iter_fake_random_record(
    n_seqs: int, size_int: Tuple[int, int], cyclic: bool
) -> List[SeqRecord]:
    for i in range(n_seqs):
        length = random.randint(*size_int)

        name = "<random record {}".format(str(uuid4()))
        rec = biopython.random_record(length, name=name, auto_annotate=True)
        rec.id = rec.name
        biopython.randomly_annotate(rec, (100, 1000))

        if cyclic:
            make_circular([rec])
        else:
            make_linear([rec])
        yield rec


def generate_fake_library(
    n_cyclic_seqs: int,
    n_linear_seqs: int,
    n_primers: int,
    cyclic_size_int: Tuple[int, int] = (3000, 10000),
    linear_size_int: Tuple[int, int] = (100, 4000),
    primer_size_int: Tuple[int, int] = (15, 60),
    return_all: bool = False,
) -> Union[List[SeqRecord], Dict[str, List[SeqRecord]]]:
    cyclic_seqs = list(
        iter_fake_random_record(n_cyclic_seqs, size_int=cyclic_size_int, cyclic=True)
    )
    linear_seqs = list(
        iter_fake_random_record(n_linear_seqs, size_int=linear_size_int, cyclic=False)
    )
    primer_seqs = list(
        iter_fake_random_record(n_primers, size_int=primer_size_int, cyclic=False)
    )
    if return_all:
        return cyclic_seqs + linear_seqs + primer_seqs
    else:
        return {"cyclic": cyclic_seqs, "linear": linear_seqs, "short": primer_seqs}


def _add_shared_sequence(
    sequences: List[SeqRecord], length: int, name="shared sequence"
):
    new_designs = []
    shared_chunk = biopython.random_record(length, name=name)
    for design in sequences:
        i = random.randint(0, len(design))
        new_designs.append(design[:i] + shared_chunk + design[i:])
    return new_designs


# TODO: replace with with MC method
def generate_fake_designs(
    n_designs: int,
    circular: bool,
    n_cyclic_seqs: int,
    n_linear_seqs: int,
    n_primers: int,
    n_primers_from_templates: int,
    design_sequence_similarity_length: int = 0,
    cyclic_size_int: Tuple[int, int] = (3000, 10000),
    linear_size_int: Tuple[int, int] = (100, 4000),
    primer_size_int: Tuple[int, int] = (15, 60),
    plasmid_size_interval: Tuple[int, int] = (5000, 10000),
    chunk_size_interval: Tuple[int, int] = (100, 3000),
    random_chunk_prob_int: Tuple[float, float] = (0, 0.5),
    random_chunk_size_int: Tuple[int, int] = (100, 1000),
):
    library_dict = generate_fake_library(
        n_cyclic_seqs=n_cyclic_seqs,
        n_linear_seqs=n_linear_seqs,
        n_primers=n_primers,
        cyclic_size_int=cyclic_size_int,
        linear_size_int=linear_size_int,
        primer_size_int=primer_size_int,
    )
    linear_seqs = library_dict["linear"]
    cyclic_seqs = library_dict["cyclic"]
    templates = cyclic_seqs + linear_seqs
    short_seqs = library_dict["short"]

    for i in range(n_primers_from_templates):
        primer = biopython.random_record_from_library(
            templates,
            circular=False,
            size_interval=(15, 100),
            max_chunks=1,
            chunk_size_interval=(15, 60),
            random_chunk_prob_int=(0, 0),
            random_chunk_size_int=(0, 0),
        )
        short_seqs.append(primer)

    # generate designs from templates or random sequence
    designs = []
    for i in range(n_designs):
        rec = biopython.random_record_from_library(
            templates,
            circular=circular,
            size_interval=plasmid_size_interval,
            chunk_size_interval=chunk_size_interval,
            random_chunk_prob_int=random_chunk_prob_int,
            random_chunk_size_int=random_chunk_size_int,
        )
        designs.append(rec)

    if design_sequence_similarity_length:
        designs = _add_shared_sequence(designs, design_sequence_similarity_length)

    if circular:
        make_circular(designs)
    else:
        make_linear(designs)

    return {
        "design": designs,
        "cyclic": cyclic_seqs,
        "linear": linear_seqs,
        "short": short_seqs,
    }
