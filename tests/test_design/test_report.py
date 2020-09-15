"""test_report.py.

Tests that generate an alignment report
"""
import pylab as plt
import pytest

from dasi import Design
from dasi.design.report import Report


class TestReport:
    @pytest.fixture(scope="module")
    def design(self):
        d = Design.fake(n_designs=3)
        d.run()
        return d

    def test_plot_coverage(self, design):
        report = Report(design)
        for qk in design.containers:
            report.plot_coverage_of_container(design.containers[qk], design.seqdb[qk])
            plt.show()

    def test_report_from_design(self, design):
        assert design.report()

    def test_report_from_design_plot_coverage(self, design):
        coverage_plots = design.report().plot_coverage()
        for fig, axes in coverage_plots.values():
            fig.show()
