import { describe, expect, it, vi } from 'vitest';
import { api } from './api';

describe('typed API client', () => {
  it('encodes variant identifiers', async () => {
    const fetchMock = vi.spyOn(globalThis, 'fetch').mockResolvedValue(
      new Response(JSON.stringify({ canonical_protein: 'F508del' }), { status: 200 })
    );
    await api.getVariant('F 508/del');
    expect(fetchMock.mock.calls[0][0]).toContain('/variants/F%20508%2Fdel');
    fetchMock.mockRestore();
  });

  it('surfaces backend errors', async () => {
    vi.spyOn(globalThis, 'fetch').mockResolvedValue(
      new Response(JSON.stringify({ detail: 'Variant not found' }), { status: 404 })
    );
    await expect(api.getVariant('nope')).rejects.toThrow('Variant not found');
    vi.restoreAllMocks();
  });

  it('fetches database statistics from /stats', async () => {
    const mockStats = {
      database_connected: true,
      evidence_loaded: true,
      total_evidence_records: 4756,
      grouped_protein_changes: 3353,
      cftr2_records: 2092,
      clinvar_records: 2664,
      cftr_structures_count: 3,
      control_structures_count: 2,
      comparable_positions: 1126,
    };
    const fetchMock = vi.spyOn(globalThis, 'fetch').mockResolvedValue(
      new Response(JSON.stringify(mockStats), { status: 200 })
    );
    const res = await api.stats();
    expect(res.total_evidence_records).toBe(4756);
    expect(res.grouped_protein_changes).toBe(3353);
    expect(fetchMock.mock.calls[0][0]).toContain('/stats');
    fetchMock.mockRestore();
  });

  it('serializes listVariants query parameters correctly', async () => {
    const fetchMock = vi.spyOn(globalThis, 'fetch').mockResolvedValue(
      new Response(JSON.stringify({ total: 1, limit: 10, offset: 0, results: [] }), { status: 200 })
    );
    await api.listVariants({ search: 'G551D', domain: 'NBD1', limit: 10 });
    const url = String(fetchMock.mock.calls[0][0]);
    expect(url).toContain('search=G551D');
    expect(url).toContain('domain=NBD1');
    expect(url).toContain('limit=10');
    fetchMock.mockRestore();
  });

  it('submits batch variants payload correctly', async () => {
    const mockBatchResp = {
      batch_id: 'batch-test-123',
      timestamp: '2026-10-12T00:00:00Z',
      total_input_rows: 2,
      unique_identifiers: 2,
      resolved_count: 2,
      unresolved_count: 0,
      evidence_available_count: 2,
      results: [
        { row_number: 1, original_input: 'G551D', status: 'RESOLVED', has_structural_evidence: true },
        { row_number: 2, original_input: 'F508del', status: 'RESOLVED', has_structural_evidence: true },
      ],
    };
    const fetchMock = vi.spyOn(globalThis, 'fetch').mockResolvedValue(
      new Response(JSON.stringify(mockBatchResp), { status: 200 })
    );
    const res = await api.batchVariants(['G551D', 'F508del']);
    expect(res.batch_id).toBe('batch-test-123');
    expect(res.resolved_count).toBe(2);
    expect(fetchMock.mock.calls[0][0]).toContain('/variants/batch');
    const init = fetchMock.mock.calls[0][1] as RequestInit;
    expect(init.method).toBe('POST');
    expect(JSON.parse(String(init.body))).toEqual({ identifiers: ['G551D', 'F508del'] });
    fetchMock.mockRestore();
  });

  it('generates structure coordinate url and residue mapping request', async () => {
    expect(api.getStructureCoordinatesUrl('5UAK')).toContain('/structures/5UAK/coordinates');

    const fetchMock = vi.spyOn(globalThis, 'fetch').mockResolvedValue(
      new Response(
        JSON.stringify({
          structure_id: '5UAK',
          canonical_position: 551,
          auth_asym_id: 'A',
          auth_seq_id: '551',
          mapping_status: 'MAPPED_OBSERVED',
        }),
        { status: 200 }
      )
    );
    const mapping = await api.getStructureResidue('5UAK', 551);
    expect(mapping.auth_asym_id).toBe('A');
    expect(mapping.mapping_status).toBe('MAPPED_OBSERVED');
    expect(fetchMock.mock.calls[0][0]).toContain('/structures/5UAK/residues/551');
    fetchMock.mockRestore();
  });
});

