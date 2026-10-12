export type VariantRecord = {
  source_db: string;
  source_record_id: string;
  source_variant_name: string | null;
  canonical_protein: string | null;
  canonical_position: number | null;
  variant_type: string;
  source_determination: string | null;
  clinvar_review_status: string | null;
  eligibility_status: string | null;
  domain_name: string | null;
  observation_pattern: string | null;
  is_fully_comparable: boolean | null;
};

export type VariantList = {
  total: number;
  limit: number;
  offset: number;
  results: VariantRecord[];
};

export type StructureEvidence = {
  state: string;
  mapping_status: string | null;
  structure_aa: string | null;
  feature_status: string | null;
  failure_reason: string | null;
  sasa_angstrom2: number | null;
  rsa_ratio: number | null;
  exposure_category: string | null;
  total_contacts: number | null;
  intra_contacts: number | null;
  inter_contacts: number | null;
  secondary_structure: string | null;
  construct_note: string | null;
  contact_definition?: string | null;
};

export type CFTR2Record = {
  legacy_name: string;
  determination_verbatim: string | null;
  alleles_reported: number | null;
  validation_flag: string;
  eligibility_status: string;
  ineligibility_reason: string | null;
};

export type ClinVarRecord = {
  variation_id: string;
  vcv_accession: string | null;
  name: string | null;
  classification_verbatim: string | null;
  review_status: string | null;
  validation_flag: string;
  eligibility_status: string;
  ineligibility_reason: string | null;
};

export type ComparativeMetrics = {
  observation_pattern: string | null;
  construct_variation: string | null;
  is_fully_comparable: boolean;
  all_structures_observed?: boolean;
  canonical_identity_match?: boolean;
  construct_context_caveat?: boolean;
  environment_comparability?: string;
  comparison_caveats?: string[];
  delta_rsa: number | null;
  delta_contacts: number | null;
  sec_struct_consistency: string | null;
  comparison_notes: string | null;
};

export type DatabaseStats = {
  database_connected: boolean;
  evidence_loaded: boolean;
  total_evidence_records: number;
  grouped_protein_changes: number;
  cftr2_records: number;
  clinvar_records: number;
  cftr_structures_count: number;
  control_structures_count: number;
  comparable_positions: number;
  error?: string;
};

export type VariantDetail = {
  canonical_protein: string | null;
  canonical_position: number | null;
  ref_aa: string | null;
  alt_aa: string | null;
  domain_name: string | null;
  source_associations: {
    cftr2_count: number;
    cftr2_records: CFTR2Record[];
    clinvar_count: number;
    clinvar_records: ClinVarRecord[];
  };
  structures: Record<'5UAK' | '6MSM' | '8EIQ', StructureEvidence>;
  comparative_metrics: ComparativeMetrics;
  warnings: { type: string; message: string }[];
  guardrail_notice: string;
};

export type BatchRowResult = {
  row_number: number;
  original_input: string;
  normalized_identifier: string | null;
  status: 'RESOLVED' | 'NOT_FOUND' | 'UNSUPPORTED_FORMAT' | 'EMPTY';
  canonical_protein: string | null;
  canonical_position: number | null;
  ref_aa: string | null;
  alt_aa: string | null;
  domain_name: string | null;
  source_record_count: number;
  source_determinations: string[];
  has_structural_evidence: boolean;
  structures: Record<'5UAK' | '6MSM' | '8EIQ', StructureEvidence> | null;
  comparative_metrics: ComparativeMetrics | null;
  warnings: { type: string; message: string }[];
  failure_reason: string | null;
};

export type BatchVariantResponse = {
  batch_id: string;
  timestamp: string;
  total_input_rows: number;
  unique_identifiers: number;
  resolved_count: number;
  unresolved_count: number;
  evidence_available_count: number;
  results: BatchRowResult[];
};

export type StructureResidueMapping = {
  structure_id: string;
  canonical_position: number;
  canonical_amino_acid: string;
  entity_id: string;
  label_asym_id: string;
  auth_asym_id: string;
  label_seq_id: number | null;
  auth_seq_id: string | null;
  insertion_code: string | null;
  structure_amino_acid: string | null;
  has_observed_residue: boolean;
  has_ca_atom: boolean;
  reference_identity_match: boolean | null;
  mapping_status: string;
  mapping_warning: string | null;
};

export type StructureMetadata = {
  id: string;
  title: string | null;
  experimental_method: string | null;
  resolution: number | null;
  structure_type: string;
  has_cftr: boolean;
  deposition_date: string | null;
  chains: {
    label_asym_id: string;
    auth_asym_id: string;
    is_cftr: boolean;
    entity_description: string | null;
    polymer_seq_len: number;
    canonical_coverage: number | null;
    observed_coverage: number | null;
  }[];
};

const base = (import.meta.env.VITE_API_BASE_URL as string | undefined) ?? '/api';

async function request<T>(path: string, init?: RequestInit): Promise<T> {
  const response = await fetch(`${base}${path}`, init);
  if (!response.ok) {
    let message = `API request failed (${response.status})`;
    try {
      const body = await response.json();
      message = body.detail ?? message;
    } catch {
      // Use fallback status message if JSON parse fails
    }
    throw new Error(message);
  }
  return response.json() as Promise<T>;
}

export const api = {
  listVariants: (params: Record<string, string | number | undefined> = {}) => {
    const query = new URLSearchParams();
    Object.entries(params).forEach(([key, value]) => {
      if (value !== undefined && value !== '') {
        query.set(key, String(value));
      }
    });
    return request<VariantList>(`/variants?${query.toString()}`);
  },
  getVariant: (identifier: string) =>
    request<VariantDetail>(`/variants/${encodeURIComponent(identifier)}`),
  batchVariants: (identifiers: string[]) =>
    request<BatchVariantResponse>('/variants/batch', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ identifiers }),
    }),
  getStructureCoordinatesUrl: (structureId: string) =>
    `${base}/structures/${encodeURIComponent(structureId)}/coordinates`,
  getStructureResidue: (structureId: string, canonicalPosition: number) =>
    request<StructureResidueMapping>(
      `/structures/${encodeURIComponent(structureId)}/residues/${canonicalPosition}`
    ),
  getStructureMetadata: (structureId: string) =>
    request<StructureMetadata>(
      `/structures/${encodeURIComponent(structureId)}/metadata`
    ),
  stats: () => request<DatabaseStats>('/stats'),
  health: () => request<{ status: string }>('/health'),
};

