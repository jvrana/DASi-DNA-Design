from dasi.alignments import Alignment
from dasi.alignments import PCRProductAlignmentGroup
from dasi.utils import Region


class TestPCRProduct:
    def test_basic_overhang(self):
        c = 3000
        fwd_query = Region(1000 - 40, 1000 + 20, c)
        fwd_subject = Region(0, 60, 60)

        rev_query = Region(2000 - 20, 2000 + 40, c)
        rev_subject = Region(0, 60, 60)

        template_query = Region(1000, 2000, c)
        template_subject = Region(0, 1000, 6000)

        fwd = Alignment(fwd_query, fwd_subject, None, None, None)
        rev = Alignment(rev_query, rev_subject, None, None, None)
        template = Alignment(template_query, template_subject, None, None, None)

        group = PCRProductAlignmentGroup(
            fwd=fwd, template=template, rev=rev, group_type=None
        )

        assert group.query_region.a == 1000 - 40
        assert group.query_region.b == 2000 + 40

    def test_primers_inside_template(self):
        c = 3000
        fwd_query = Region(1000 - 40, 1000 + 20, c)
        fwd_subject = Region(0, 60, 60)

        rev_query = Region(2000 - 20, 2000 + 40, c)
        rev_subject = Region(0, 60, 60)

        template_query = Region(500, 2500, c)
        template_subject = Region(0, 2000, 6000)

        fwd = Alignment(fwd_query, fwd_subject, None, None, None)
        rev = Alignment(rev_query, rev_subject, None, None, None)
        template = Alignment(template_query, template_subject, None, None, None)

        group = PCRProductAlignmentGroup(
            fwd=fwd, template=template, rev=rev, group_type=None
        )

        assert group.query_region.a == 1000 - 40
        assert group.query_region.b == 2000 + 40


# class TestMultiPCRPRoduct:
#
#     def test_basic_get_template
