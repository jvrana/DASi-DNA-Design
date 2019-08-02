
class AssemblyGraphBuilder(object):
    def __init__(self, alignment_container: AlignmentContainer):
        self.container = alignment_container
        self.span_cost = SpanCost()
        self.G = None

    def _edge_weight(self, type1: str, type2: str, span: int) -> float:
        ext = self._extension_tuple_from_types(type1, type2)
        assembly_cost = np.clip(self.span_cost.cost(span, ext), 0, Constants.INF)
        frag_cost1 = Constants.PCR_COST[type1] # we only count the source fragment cost
        return assembly_cost + frag_cost1

    def _extension_tuple_from_types(self, left_type, right_type):
        """
        Converts two "TYPES" into an 'extension tuple' that can be consumed by
        the cost classes. The 'extension tuple' is a tuple representing the ability
        of the fragments across a junction to be design to bridge a gap. For example:

        (0, 0) means the following gap cannot be spanned:

        ::
            (0, 0)
            ------|      |------

        It is also important to note that "0" means the either the primer for that fragment
        exists

        While the follow gap *may* be spanned by designing the left fragment with a longer
        primer:

        ::

            (1, 0)
            -------|    |-------
               <----------

        And finally, the following gap may be spanned by both primers.


        ::

            (1, 1)
                    -------->
            -------|    |-------
               <----------

        :param left_type:
        :param right_type:
        :return:
        """
        if (
            left_type == Constants.PCR_PRODUCT_WITH_PRIMERS
            or Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER
        ):
            left_extendable = 1
        else:
            left_extendable = 0

        if (
            right_type == Constants.PCR_PRODUCT_WITH_PRIMERS
            or Constants.PCR_PRODUCT_WITH_LEFT_PRIMER
        ):
            right_extendable = 1
        else:
            right_extendable = 0
        return (left_extendable, right_extendable)

    def _add_node(self, g1: AlignmentGroup):
        """Add a node from an alignment group to the graph."""
        n1 = self.container.alignment_hash(g1.alignments[0])
        self.G.add_node(
            n1,
            start=g1.query_region.start,
            end=g1.query_region.end,
            region=str(g1.query_region),
        )
        return n1

    def _add_edge(self, g1, g2, span_length, **kwargs):
        """Add nodes, edge, and edge_weight from two alignment groups. A spanning distance (bp)
        must be provided."""
        n1 = self._add_node(g1)
        n2 = self._add_node(g2)
        edata = {
            "weight": self._edge_weight(g1.type, g2.type, span_length),
            "region_1": str(g1.query_region),
            "region_2": str(g2.query_region),
        }
        edata.update(kwargs)
        self.G.add_edge(n1, n2, **edata)

    def add_edge(self, group, other_group):
        """
        Create an edge between one alignment group and another.

        :param group:
        :param other_group:
        :return:
        """
        # verify contexts are the same
        r1 = group.query_region
        r2 = other_group.query_region
        assert r2.same_context(r1)

        if r1 in r2 or r2 in r1:
            return
        if r2.a in r1 and r2.b not in r1:
            # strict forward overlap
            overlap = r1.intersection(r2)
            if not overlap:
                raise Exception("We expected an overlap here")
            # TODO: penalize small overlap
            self._add_edge(
                group, other_group, span_length=-len(overlap), name="overlap"
            )
        elif r2.a not in r1 and r2.b not in r1:
            try:
                connecting_span = r1.connecting_span(r2)
            except Exception as e:
                raise e
            if connecting_span:
                span_length = len(connecting_span)
            elif r1.consecutive(r2):
                span_length = 0
            else:
                return
            self._add_edge(
                group, other_group, span_length=span_length, name="synthesis"
            )

    # TODO: linear assemblies
    # TODO: replace region methods with something faster.
    # TODO: edge cost should include (i) cost of producing the source and (ii) cost of assembly
    # Verify queries have same context
    def build_assembly_graph(self):
        self.G = nx.DiGraph(name="Assembly Graph")
        # TODO: tests for validating over-origin edges are being produced
        # produce non-spanning edges
        # makes edges for any regions
        groups, group_keys = sort_with_keys(
            self.container.get_groups_by_types(
                [
                    Constants.PCR_PRODUCT,
                    Constants.PCR_PRODUCT_WITH_RIGHT_PRIMER,
                    Constants.PCR_PRODUCT_WITH_LEFT_PRIMER,
                    Constants.PCR_PRODUCT_WITH_PRIMERS,
                    Constants.FRAGMENT,
                ]
            ),
            key=lambda g: g.query_region.a,
        )
        for g1, g2 in itertools.product(groups, repeat=2):
            self.add_edge(g1, g2)
        return self.G
