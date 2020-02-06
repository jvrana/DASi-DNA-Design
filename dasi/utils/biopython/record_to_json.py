from copy import deepcopy

from Bio.SeqFeature import CompoundLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord


def feature_to_json(feature: SeqFeature, length: int) -> dict:
    location = feature.location
    if isinstance(location, CompoundLocation):
        start = location.parts[0].start
        end = location.parts[-1].end
        strand = location.parts[0].strand
        assert location.parts[0].end == length
    else:
        start = int(feature.location.start)
        end = int(feature.location.end)
        strand = feature.location.strand
    qualifier_label = feature.qualifiers["label"]
    if isinstance(qualifier_label, list):
        name = qualifier_label[0]
    elif isinstance(qualifier_label, str):
        name = qualifier_label
    elif qualifier_label is None:
        name = "misc_feature"
    else:
        raise ValueError(
            "Type '{}' {} is not a valid feature name".format(
                type(qualifier_label), qualifier_label
            )
        )
    color = feature.qualifiers["ApEinfo_fwdcolor"][0]
    if end >= length:
        end = 0
    if start >= length:
        start = 0
    return {
        "start": start,
        "end": end,
        "strand": strand,
        "color": color,
        "name": name,
        "type": feature.type,
    }


def seqrecord_to_json(record: SeqRecord) -> dict:
    feature_annotations = [feature_to_json(f, len(record.seq)) for f in record.features]
    data = {
        "bases": str(record.seq),
        "length": len(record.seq),
        "name": record.name,
        "id": record.id,
        "annotations": feature_annotations,
        "customFields": deepcopy(record.annotations),
    }
    annotations = deepcopy(record.annotations)

    if annotations.get("topology", "linear") == "circular":
        annotations["isCircular"] = True
    else:
        annotations["isCircular"] = False
    if "topology" in annotations:
        del annotations["topology"]
    data.update(annotations)
    return data
