"""Tests generation of fake SeqRecords."""
from dasi.utils.sequence.sequence_generator import generate_fake_designs


def test_fake_designs():
    print(generate_fake_designs(3, True, 50, 50, 50, 50))


def test_fake_similar_designs():
    print(
        generate_fake_designs(
            3, True, 50, 50, 50, 50, design_sequence_similarity_length=1000
        )
    )
