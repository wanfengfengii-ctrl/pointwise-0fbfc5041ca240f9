from .result_store import ResultStore
from .result_codec import SUPERCELL_CALCULATION, result_to_mapping
from .supercell_codec import supercell_to_mapping
from .export import (
    build_export_packet, build_supercell_export_packet, encode_export_packet)


class Results:
    def __init__(self, archive):
        self.store = ResultStore(archive.connection)
        self.revisions = archive.revisions

    def read(self, result_id):
        result = self.store.read_bundle(result_id)
        if result.source[3] == SUPERCELL_CALCULATION:
            return supercell_to_mapping(result)
        return result_to_mapping(result)

    def export(self, result_id):
        result = self.store.read_bundle(result_id)
        revision = self.revisions.read_revision(
            result.source[0], result.source[1])
        if result.source[3] == SUPERCELL_CALCULATION:
            packet = build_supercell_export_packet(revision, result)
        else:
            packet = build_export_packet(revision, result)
        return encode_export_packet(packet)
