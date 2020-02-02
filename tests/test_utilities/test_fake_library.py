"""Tests generation of fake SeqRecords."""
from dasi.utils.sequence_generator import fake_designs


def test_fake_designs():
    print(fake_designs(3, True, 50, 50, 50, 50))


def test_fake_similar_designs():
    print(fake_designs(3, True, 50, 50, 50, 50, design_sequence_similarity_length=1000))
