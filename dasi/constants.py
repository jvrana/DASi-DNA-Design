
class Constants(object):
    FRAGMENT = (
        "PRE-MADE DNA FRAGMENT"
    )  # an alignment that is generate from an already existing PCR product or fragment
    PCR_PRODUCT = (
        "PCR_PRODUCT"
    )  # an alignment that is to be generated from a PCR product
    PCR_PRODUCT_WITH_PRIMERS = (
        "PCR_PRODUCT_WITH_PRIMERS"
    )  # PCR product that can be produces from existing primers
    PCR_PRODUCT_WITH_LEFT_PRIMER = (
        "PCR_PRODUCT_WITH_LEFT_PRIMER"
    )  # PCR product with existing left primer
    PCR_PRODUCT_WITH_RIGHT_PRIMER = (
        "PCR_PRODUCT_WITH_RIGHT_PRIMER"
    )  # PCR product with existing right primer
    SHARED_FRAGMENT = (
        "FRAGMENT_SHARED_WITH_OTHER_QUERIES"
    )  # A fragment alignment that is shared with other queries for potential reuse

    PCR_COST = {
        PCR_PRODUCT: 30 + 25 + 10,
        PCR_PRODUCT_WITH_PRIMERS: 30,
        PCR_PRODUCT_WITH_LEFT_PRIMER: 15,
        PCR_PRODUCT_WITH_RIGHT_PRIMER: 15,
        FRAGMENT: 0,
    }

    PRIMER = "PRIMER"  # a primer binding alignment

    PRIMER_MIN_BIND = 15
    MIN_OVERLAP = 20
    MAX_HOMOLOGY = 100
    INF = 10.0 ** 6

