from .exact_numbers import json_dump, json_load
from .input_model import cell_to_mapping
from .input_digest import verify_cell_digest
from .intake import accept_cell
from .result_codec import encode_result
from .executor import validate_geometry_result
from .supercell_codec import encode_supercell_result
from .supercell_audit import audit_supercell_evidence

EXPORT_FORMAT_VERSION = 1
SUPERCELL_EXPORT_FORMAT = "crystal-supercell-evidence"


def validate_export_evidence(revision, result):
    if result.source[:3] != (
            revision.product_id, revision.revision, revision.digest):
        raise ValueError("export source identity mismatch")
    verify_cell_digest(revision.cell, revision.digest)
    raw_cell, raw_digest = accept_cell(revision.cell.raw_mapping())
    if raw_digest != revision.digest:
        raise ValueError("raw and canonical evidence differ")
    validate_geometry_result(raw_cell, result)
    return result


def validate_supercell_export_evidence(revision, result):
    if result.source[:3] != (
            revision.product_id, revision.revision, revision.digest):
        raise ValueError("export source identity mismatch")
    verify_cell_digest(revision.cell, revision.digest)
    raw_cell, raw_digest = accept_cell(revision.cell.raw_mapping())
    if raw_digest != revision.digest:
        raise ValueError("raw and canonical evidence differ")
    audit_supercell_evidence(raw_cell, result)
    return result


def _source_mapping(revision):
    return {
        "product_id": revision.product_id, "revision": revision.revision,
        "digest": revision.digest,
        "cell": cell_to_mapping(revision.cell),
        "raw_cell": revision.cell.raw_mapping(),
    }


def build_export_packet(revision, result):
    validate_export_evidence(revision, result)
    return {
        "format": "crystal-primitive-evidence",
        "version": EXPORT_FORMAT_VERSION,
        "source": _source_mapping(revision),
        "result": json_load(encode_result(result))["result"],
    }


def build_supercell_export_packet(revision, result):
    validate_supercell_export_evidence(revision, result)
    return {
        "format": SUPERCELL_EXPORT_FORMAT,
        "version": EXPORT_FORMAT_VERSION,
        "source": _source_mapping(revision),
        "result": json_load(encode_supercell_result(result))["result"],
    }


def encode_export_packet(packet):
    return json_dump(packet)
