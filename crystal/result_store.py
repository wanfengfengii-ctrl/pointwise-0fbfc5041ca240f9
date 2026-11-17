from .database import transaction
from .exact_numbers import read_integer
from .revision_store import RevisionStore
from .result_codec import (
    CALCULATION, SUPERCELL_CALCULATION, encode_result, decode_result,
    source_to_mapping)
from .supercell_codec import encode_supercell_result, decode_supercell_result
from .supercell_audit import audit_supercell_evidence


ENCODERS = {
    CALCULATION: encode_result,
    SUPERCELL_CALCULATION: encode_supercell_result,
}
DECODERS = {
    CALCULATION: decode_result,
    SUPERCELL_CALCULATION: decode_supercell_result,
}


class ResultStore:
    def __init__(self, connection):
        self.connection = connection
        self.revisions = RevisionStore(connection)

    def _publish(self, task, result, calculation, encode):
        source = task.source
        source_to_mapping(source)
        if source[3] != calculation:
            raise ValueError("task calculation differs from its publisher")
        if result.source != source:
            raise ValueError("task/result source mismatch")
        payload = encode(result)
        with transaction(self.connection):
            revision = self.revisions.read_revision(source[0], source[1])
            if revision.digest != source[2]:
                raise ValueError("task/revision source mismatch")
            cursor = self.connection.execute(
                "INSERT INTO tasks(product_id,revision,digest,calculation) "
                "VALUES (?,?,?,?)", source)
            task_id = cursor.lastrowid
            cursor = self.connection.execute(
                "INSERT INTO results(task_id) VALUES (?)", (task_id,))
            result_id = cursor.lastrowid
            self.connection.execute(
                "INSERT INTO result_bodies(result_id,payload) VALUES (?,?)",
                (result_id, payload))
        return result_id

    def publish_bundle(self, task, result):
        return self._publish(task, result, CALCULATION, encode_result)

    def publish_supercell_bundle(self, task, result):
        return self._publish(task, result, SUPERCELL_CALCULATION,
                             encode_supercell_result)

    def read_bundle(self, result_id):
        read_integer(result_id, "result_id")
        if result_id < 1:
            raise ValueError("result_id must be positive")
        row = self.connection.execute(
            "SELECT t.product_id,t.revision,t.digest,t.calculation,b.payload "
            "FROM results r JOIN tasks t ON t.task_id=r.task_id "
            "LEFT JOIN result_bodies b ON b.result_id=r.result_id "
            "WHERE r.result_id=?", (result_id,)).fetchone()
        if row is None:
            raise KeyError(result_id)
        if row[4] is None:
            raise ValueError("result body is missing")
        calculation = row[3]
        if calculation not in DECODERS:
            raise ValueError("unsupported stored calculation")
        result = DECODERS[calculation](row[4])
        if result.source != tuple(row[:4]):
            raise ValueError("result header/body binding mismatch")
        revision = self.revisions.read_revision(row[0], row[1])
        if revision.digest != row[2]:
            raise ValueError("result/revision binding mismatch")
        if calculation == SUPERCELL_CALCULATION:
            audit_supercell_evidence(revision.cell, result)
        return result
