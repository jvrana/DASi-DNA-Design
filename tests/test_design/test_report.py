import pylab as plt

from dasi.design.report import Report


class TestReport:
    def test_plot_coverage(self, single_processed_results):
        design, result = single_processed_results
        report = Report(design)
        for qk in design.containers:
            report.plot_coverage(design.containers[qk], design.seqdb[qk])
            plt.show()
