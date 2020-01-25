"""Tests generation of fake SeqRecords."""
from dasi.utils.testing_utils import fake_designs


def test_fake_designs():
    print(fake_designs(3, True, 50, 50, 50, 50))
