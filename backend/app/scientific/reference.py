"""Canonical CFTR sequence reference and residue validation.

Reference sequence:
    UniProtKB: P13569 (CFTR_HUMAN)
    Entry version: SV=3
    Sequence length: 1480 amino acids
"""

from dataclasses import dataclass
from typing import Optional

UNIPROT_ACCESSION = "P13569"
UNIPROT_ID = "CFTR_HUMAN"
SEQUENCE_LENGTH = 1480
SEQUENCE_SHA256 = "a86f7dbeb9c5f48401a1a452b440d782c35067d2da152180a09ce26f76d11639"

CFTR_CANONICAL_SEQUENCE = (
    "MQRSPLEKASVVSKLFFSWTRPILRKGYRQRLELSDIYQIPSVDSADNLSEKLEREWDRE"
    "LASKKNPKLINALRRCFFWRFMFYGIFLYLGEVTKAVQPLLLGRIIASYDPDNKEERSIA"
    "IYLGIGLCLLFIVRTLLLHPAIFGLHHIGMQMRIAMFSLIYKKTLKLSSRVLDKISIGQL"
    "VSLLSNNLNKFDEGLALAHFVWIAPLQVALLMGLIWELLQASAFCGLGFLIVLALFQAGL"
    "GRMMMKYRDQRAGKISERLVITSEMIENIQSVKAYCWEEAMEKMIENLRQTELKLTRKAA"
    "YVRYFNSSAFFFSGFFVVFLSVLPYALIKGIILRKIFTTISFCIVLRMAVTRQFPWAVQT"
    "WYDSLGAINKIQDFLQKQEYKTLEYNLTTTEVVMENVTAFWEEGFGELFEKAKQNNNNRK"
    "TSNGDDSLFFSNFSLLGTPVLKDINFKIERGQLLAVAGSTGAGKTSLLMVIMGELEPSEG"
    "KIKHSGRISFCSQFSWIMPGTIKENIIFGVSYDEYRYRSVIKACQLEEDISKFAEKDNIV"
    "LGEGGITLSGGQRARISLARAVYKDADLYLLDSPFGYLDVLTEKEIFESCVCKLMANKTR"
    "ILVTSKMEHLKKADKILILHEGSSYFYGTFSELQNLQPDFSSKLMGCDSFDQFSAERRNS"
    "ILTETLHRFSLEGDAPVSWTETKKQSFKQTGEFGEKRKNSILNPINSIRKFSIVQKTPLQ"
    "MNGIEEDSDEPLERRLSLVPDSEQGEAILPRISVISTGPTLQARRRQSVLNLMTHSVNQG"
    "QNIHRKTTASTRKVSLAPQANLTELDIYSRRLSQETGLEISEEINEEDLKECFFDDMESI"
    "PAVTTWNTYLRYITVHKSLIFVLIWCLVIFLAEVAASLVVLWLLGNTPLQDKGNSTHSRN"
    "NSYAVIITSTSSYYVFYIYVGVADTLLAMGFFRGLPLVHTLITVSKILHHKMLHSVLQAP"
    "MSTLNTLKAGGILNRFSKDIAILDDLLPLTIFDFIQLLLIVIGAIAVVAVLQPYIFVATV"
    "PVIVAFIMLRAYFLQTSQQLKQLESEGRSPIFTHLVTSLKGLWTLRAFGRQPYFETLFHK"
    "ALNLHTANWFLYLSTLRWFQMRIEMIFVIFFIAVTFISILTTGEGEGRVGIILTLAMNIM"
    "STLQWAVNSSIDVDSLMRSVSRVFKFIDMPTEGKPTKSTKPYKNGQLSKVMIIENSHVKK"
    "DDIWPSGGQMTVKDLTAKYTEGGNAILENISFSISPGQRVGLLGRTGSGKSTLLSAFLRL"
    "LNTEGEIQIDGVSWDSITLQQWRKAFGVIPQKVFIFSGTFRKNLDPYEQWSDQEIWKVAD"
    "EVGLRSVIEQFPGKLDFVLVDGGCVLSHGHKQLMCLARSVLSKAKILLLDEPSAHLDPVT"
    "YQIIRRTLKQAFADCTVILCEHRIEAMLECQQFLVIEENKVRQYDSIQKLLNERSLFRQA"
    "ISPSDRVKLFPHRNSSKCKSKPQIAALKEETEEEVQDTRL"
)

# Standard 20 amino acid codes mapping (3-letter to 1-letter)
AA3_TO_AA1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O",
}

# 1-letter to 3-letter
AA1_TO_AA3 = {v: k for k, v in AA3_TO_AA1.items()}


@dataclass(frozen=True)
class ValidationResult:
    """Result of validating a reported residue against canonical UniProt P13569."""
    is_valid: bool
    pos: int
    reported_ref: str
    canonical_ref: Optional[str]
    flag: str
    message: str


def get_canonical_residue(pos: int) -> Optional[str]:
    """Retrieve 1-based canonical amino acid from UniProt P13569 sequence.
    
    Args:
        pos: 1-based amino acid position.
        
    Returns:
        Single-letter amino acid code, or None if pos is out of bounds (1..1480).
    """
    if 1 <= pos <= SEQUENCE_LENGTH:
        return CFTR_CANONICAL_SEQUENCE[pos - 1]
    return None


def validate_wt_residue(pos: int, reported_ref: str) -> ValidationResult:
    """Validate that reported reference amino acid matches UniProt P13569 at pos.
    
    CRITICAL SCIENTIFIC CONSTRAINT:
    Mismatches must be explicitly flagged and NEVER silently corrected.
    
    Args:
        pos: 1-based amino acid position in CFTR sequence.
        reported_ref: Single-letter (or 3-letter) reported reference amino acid.
        
    Returns:
        ValidationResult indicating match, mismatch, or out of bounds.
    """
    # Normalize 3-letter code if provided
    ref = reported_ref.strip()
    if len(ref) == 3:
        ref = AA3_TO_AA1.get(ref.capitalize(), ref)
    ref = ref.upper()

    if pos < 1 or pos > SEQUENCE_LENGTH:
        return ValidationResult(
            is_valid=False,
            pos=pos,
            reported_ref=reported_ref,
            canonical_ref=None,
            flag="OUT_OF_BOUNDS",
            message=f"Position {pos} is outside canonical CFTR range (1..{SEQUENCE_LENGTH}).",
        )

    canonical = CFTR_CANONICAL_SEQUENCE[pos - 1]
    if ref == canonical:
        return ValidationResult(
            is_valid=True,
            pos=pos,
            reported_ref=reported_ref,
            canonical_ref=canonical,
            flag="MATCH",
            message=f"Residue {ref} matches canonical CFTR at position {pos}.",
        )
    else:
        # Contextual explanation for the known historical M470 polymorphism
        extra_note = ""
        if pos == 470 and ref == "M" and canonical == "V":
            extra_note = (
                " Note: Historical CFTR cDNAs (e.g. GenBank M28668) contained 1540A encoding Met470, "
                "leading to legacy nomenclature 'M470V'. Modern canonical UniProt P13569 and RefSeq NM_000492 "
                "standardize on Val470 (rs1800089). CFTR2 curation notes confirm genomic coordinates are omitted "
                "for M470V because canonical reference has V."
            )
        return ValidationResult(
            is_valid=False,
            pos=pos,
            reported_ref=reported_ref,
            canonical_ref=canonical,
            flag="WT_MISMATCH",
            message=(
                f"Reference mismatch at position {pos}: reported '{reported_ref}' ({ref}), "
                f"expected canonical '{canonical}'.{extra_note}"
            ),
        )
