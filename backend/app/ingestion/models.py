"""Ingestion data models and summary structures."""

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, Optional


@dataclass
class ImportSummary:
    """Summary of ingestion execution for reproducible audit trail."""
    source_name: str
    source_file: str
    file_sha256: str
    source_version: Optional[str]
    timestamp: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    rows_read: int = 0
    valid_rows: int = 0
    skipped_rows: int = 0
    duplicate_identifiers: int = 0
    normalization_failures: int = 0
    wt_matches: int = 0
    wt_mismatches: int = 0
    wt_out_of_bounds: int = 0
    unvalidated_rows: int = 0
    exclusion_reasons: Dict[str, int] = field(default_factory=dict)

    def add_exclusion_reason(self, reason: str, count: int = 1) -> None:
        self.exclusion_reasons[reason] = self.exclusion_reasons.get(reason, 0) + count

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def to_text_report(self) -> str:
        lines = [
            f"=== {self.source_name} Ingestion Summary ===",
            f"Source file:            {self.source_file}",
            f"SHA-256:                {self.file_sha256}",
            f"Source version:         {self.source_version or 'N/A'}",
            f"Timestamp (UTC):        {self.timestamp}",
            f"Total rows read:        {self.rows_read}",
            f"Valid rows imported:    {self.valid_rows}",
            f"Skipped rows:           {self.skipped_rows}",
            f"Duplicate identifiers:  {self.duplicate_identifiers}",
            f"Normalization failures: {self.normalization_failures}",
            f"UniProt WT matches:     {self.wt_matches}",
            f"UniProt WT mismatches:  {self.wt_mismatches}",
            f"Out of bounds pos:      {self.wt_out_of_bounds}",
            f"Unvalidated rows:       {self.unvalidated_rows}",
        ]
        if self.exclusion_reasons:
            lines.append("Exclusion breakdown:")
            for reason, count in sorted(self.exclusion_reasons.items(), key=lambda x: -x[1]):
                lines.append(f"  - {reason}: {count}")
        return "\n".join(lines)
