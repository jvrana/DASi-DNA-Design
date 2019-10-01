from pyblast import __version__
from pyblast import BioBlastFactory
from pyblast.utils import load_fasta_glob
from pyblast.utils import load_genbank_glob
from pyblast.utils import make_circular
from pyblast.utils import make_linear
from pyblast.utils import reverse_complement


PRIMERS = "primers"
TEMPLATES = "templates"
QUERIES = "queries"


def test(paths):
    factory = BioBlastFactory()

    templates = load_genbank_glob(paths[TEMPLATES])
    subject = templates[0]
    query = subject.reverse_complement()

    templates = make_linear([subject])
    queries = make_linear([query])

    factory.add_records(templates, TEMPLATES)
    factory.add_records(queries, QUERIES)

    blast = factory(TEMPLATES, QUERIES)

    results = blast.quick_blastn()
    assert results[0]["subject"]["strand"] == -1
    assert results[0]["subject"]["start"] == len(templates[0])
