import React, { useState } from 'react';
import {
  Activity,
  AlertTriangle,
  ArrowRight,
  CheckCircle2,
  ChevronLeft,
  ChevronRight,
  Download,
  FileSpreadsheet,
  FileText,
  Layers,
  Search,
  Upload,
} from 'lucide-react';
import { api, BatchRowResult, BatchVariantResponse } from '../lib/api';

const EXAMPLE_VARIANTS = [
  'G551D',
  'F508del',
  'p.Gly551Asp',
  'E1371Q',
  'M470V',
  'p.Phe508del',
  'UnknownXYZ999',
  'G551D', // duplicate test
];

function escapeCsvField(val: unknown): string {
  if (val === null || val === undefined) return '';
  const str = String(val);
  if (str.includes(',') || str.includes('"') || str.includes('\n')) {
    return `"${str.replace(/"/g, '""')}"`;
  }
  return str;
}

export function BatchAnalysis({
  onSelectVariant,
  onInspect3D,
}: {
  onSelectVariant: (variantId: string) => void;
  onInspect3D: (variantId: string, position?: number) => void;
}) {
  const [inputMode, setInputMode] = useState<'text' | 'csv'>('text');
  const [textInput, setTextInput] = useState('');
  const [csvFileName, setCsvFileName] = useState<string | null>(null);
  const [csvIdentifiers, setCsvIdentifiers] = useState<string[]>([]);
  const [csvError, setCsvError] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [response, setResponse] = useState<BatchVariantResponse | null>(null);

  // Table filtering and pagination
  const [statusFilter, setStatusFilter] = useState<'ALL' | 'RESOLVED' | 'UNRESOLVED'>('ALL');
  const [tableSearch, setTableSearch] = useState('');
  const [page, setPage] = useState(0);
  const pageSize = 10;

  // Active list of identifiers based on input mode
  const currentIdentifiers =
    inputMode === 'text'
      ? textInput
          .split('\n')
          .map((s) => s.trim())
          .filter(Boolean)
      : csvIdentifiers;

  const handleFileUpload = (e: React.ChangeEvent<HTMLInputElement>) => {
    setCsvError(null);
    const file = e.target.files?.[0];
    if (!file) return;

    if (file.size > 2 * 1024 * 1024) {
      setCsvError('File size exceeds the 2 MB limit.');
      return;
    }

    setCsvFileName(file.name);
    const reader = new FileReader();
    reader.onload = (event) => {
      try {
        const text = event.target?.result as string;
        const lines = text.split(/\r?\n/).map((l) => l.trim()).filter(Boolean);
        if (lines.length === 0) {
          setCsvError('Uploaded CSV file is empty.');
          return;
        }

        // Parse header and find variant identifier column
        const header = lines[0].split(',').map((h) => h.replace(/^["']|["']$/g, '').trim().toLowerCase());
        let colIdx = header.findIndex((h) =>
          ['variant', 'identifier', 'variant_identifier', 'protein', 'protein_change', 'name', 'mutation'].includes(h)
        );

        const startIndex = colIdx !== -1 ? 1 : 0;
        if (colIdx === -1) colIdx = 0; // Default to first column if no recognized header

        const ids: string[] = [];
        for (let i = startIndex; i < lines.length; i++) {
          const parts = lines[i].split(',').map((p) => p.replace(/^["']|["']$/g, '').trim());
          if (parts[colIdx]) {
            ids.push(parts[colIdx]);
          }
        }

        if (ids.length === 0) {
          setCsvError('Could not extract any variant identifiers from column.');
          return;
        }

        if (ids.length > 500) {
          setCsvError(`CSV contains ${ids.length} rows; batch limit is 500 rows.`);
          return;
        }

        setCsvIdentifiers(ids);
      } catch (err: unknown) {
        const msg = err instanceof Error ? err.message : 'Error parsing CSV file';
        setCsvError(msg);
      }
    };
    reader.readAsText(file);
  };

  const handleRunBatch = async () => {
    setError(null);
    const ids =
      inputMode === 'text'
        ? textInput.split(/\r?\n/)
        : csvIdentifiers;

    if (!ids || ids.length === 0 || ids.every((s) => !s.trim())) {
      setError('Please provide at least one variant identifier.');
      return;
    }

    if (ids.length > 500) {
      setError(`Batch size (${ids.length}) exceeds the maximum limit of 500 identifiers.`);
      return;
    }

    setLoading(true);
    try {
      const res = await api.batchVariants(ids);
      setResponse(res);
      setPage(0);
    } catch (err: unknown) {
      const msg = err instanceof Error ? err.message : 'Batch analysis request failed';
      setError(msg);
    } finally {
      setLoading(false);
    }
  };

  const handleLoadExamples = () => {
    setInputMode('text');
    setTextInput(EXAMPLE_VARIANTS.join('\n'));
    setError(null);
  };

  const handleDownloadTemplate = () => {
    const content = 'variant_identifier\nG551D\nF508del\np.Gly551Asp\nE1371Q\nM470V\n';
    const blob = new Blob([content], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.setAttribute('download', 'cftr_batch_template.csv');
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  };

  const handleExportResults = () => {
    if (!response || !response.results.length) return;

    const headers = [
      'Row Number',
      'Original Input',
      'Status',
      'Normalized Protein',
      'Canonical Position',
      'Domain',
      'Source Records Count',
      'Source Determinations',
      'Structural Evidence Available',
      '5UAK State',
      '5UAK SASA (A^2)',
      '5UAK RSA (%)',
      '5UAK Contacts',
      '5UAK Sec Struct',
      '6MSM State',
      '6MSM SASA (A^2)',
      '6MSM RSA (%)',
      '6MSM Contacts',
      '6MSM Sec Struct',
      '8EIQ State',
      '8EIQ SASA (A^2)',
      '8EIQ RSA (%)',
      '8EIQ Contacts',
      '8EIQ Sec Struct',
      'Delta RSA (%)',
      'Delta Contacts',
      'Sec Struct Consistency',
      'Warnings',
      'Failure Reason',
    ];

    const rows = response.results.map((r) => {
      const s5 = r.structures?.['5UAK'];
      const s6 = r.structures?.['6MSM'];
      const s8 = r.structures?.['8EIQ'];
      const m = r.comparative_metrics;

      return [
        r.row_number,
        r.original_input,
        r.status,
        r.normalized_identifier ?? '',
        r.canonical_position ?? '',
        r.domain_name ?? '',
        r.source_record_count,
        r.source_determinations.join('; '),
        r.has_structural_evidence ? 'YES' : 'NO',
        s5?.mapping_status ?? '',
        s5?.sasa_angstrom2 !== null && s5?.sasa_angstrom2 !== undefined ? s5.sasa_angstrom2.toFixed(1) : '',
        s5?.rsa_ratio !== null && s5?.rsa_ratio !== undefined ? (s5.rsa_ratio * 100).toFixed(1) : '',
        s5?.total_contacts ?? '',
        s5?.secondary_structure ?? '',
        s6?.mapping_status ?? '',
        s6?.sasa_angstrom2 !== null && s6?.sasa_angstrom2 !== undefined ? s6.sasa_angstrom2.toFixed(1) : '',
        s6?.rsa_ratio !== null && s6?.rsa_ratio !== undefined ? (s6.rsa_ratio * 100).toFixed(1) : '',
        s6?.total_contacts ?? '',
        s6?.secondary_structure ?? '',
        s8?.mapping_status ?? '',
        s8?.sasa_angstrom2 !== null && s8?.sasa_angstrom2 !== undefined ? s8.sasa_angstrom2.toFixed(1) : '',
        s8?.rsa_ratio !== null && s8?.rsa_ratio !== undefined ? (s8.rsa_ratio * 100).toFixed(1) : '',
        s8?.total_contacts ?? '',
        s8?.secondary_structure ?? '',
        m?.delta_rsa !== null && m?.delta_rsa !== undefined ? (m.delta_rsa * 100).toFixed(1) : '',
        m?.delta_contacts ?? '',
        m?.sec_struct_consistency ?? '',
        r.warnings.map((w) => `${w.type}: ${w.message}`).join('; '),
        r.failure_reason ?? '',
      ]
        .map(escapeCsvField)
        .join(',');
    });

    const csvContent = [headers.join(','), ...rows].join('\r\n');
    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.setAttribute('download', `cftr_structmap_batch_${response.batch_id}.csv`);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  };

  // Filtered rows for UI table
  const filteredResults = (response?.results || []).filter((r) => {
    if (statusFilter === 'RESOLVED' && r.status !== 'RESOLVED') return false;
    if (statusFilter === 'UNRESOLVED' && r.status === 'RESOLVED') return false;
    if (tableSearch) {
      const q = tableSearch.toLowerCase();
      const matchInput = r.original_input.toLowerCase().includes(q);
      const matchNorm = (r.normalized_identifier || '').toLowerCase().includes(q);
      const matchDomain = (r.domain_name || '').toLowerCase().includes(q);
      return matchInput || matchNorm || matchDomain;
    }
    return true;
  });

  const totalPages = Math.ceil(filteredResults.length / pageSize);
  const paginatedRows = filteredResults.slice(page * pageSize, (page + 1) * pageSize);

  return (
    <>
      <section className="page-intro">
        <span className="eyebrow">HIGH-THROUGHPUT LOOKUP</span>
        <h2>Batch Variant Analysis</h2>
        <p>
          Submit multi-variant cohorts to retrieve normalized clinical classifications,
          per-structure experimental features, and cross-structure comparative metrics.
        </p>
      </section>

      {/* Input Configuration Card */}
      <article className="card" style={{ marginBottom: '16px' }}>
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '12px' }}>
          <div style={{ display: 'flex', gap: '8px' }}>
            <button
              className={inputMode === 'text' ? 'primary' : 'secondary'}
              style={{ padding: '6px 12px', fontSize: '13px' }}
              onClick={() => setInputMode('text')}
            >
              <FileText style={{ width: '14px', marginRight: '6px', display: 'inline-block' }} />
              Paste Identifiers
            </button>
            <button
              className={inputMode === 'csv' ? 'primary' : 'secondary'}
              style={{ padding: '6px 12px', fontSize: '13px' }}
              onClick={() => setInputMode('csv')}
            >
              <FileSpreadsheet style={{ width: '14px', marginRight: '6px', display: 'inline-block' }} />
              Upload CSV File
            </button>
          </div>

          <div style={{ display: 'flex', gap: '8px' }}>
            <button
              className="secondary"
              style={{ padding: '5px 10px', fontSize: '12px' }}
              onClick={handleLoadExamples}
            >
              Load Test Examples
            </button>
            <button
              className="secondary"
              style={{ padding: '5px 10px', fontSize: '12px' }}
              onClick={handleDownloadTemplate}
            >
              <Download style={{ width: '12px', marginRight: '4px', display: 'inline-block' }} />
              Template CSV
            </button>
          </div>
        </div>

        {inputMode === 'text' ? (
          <div>
            <textarea
              rows={6}
              style={{
                width: '100%',
                padding: '10px',
                fontFamily: 'monospace',
                fontSize: '13px',
                borderRadius: '6px',
                border: '1px solid #cbd5e1',
                resize: 'vertical',
              }}
              placeholder={`Enter one variant per line (up to 500 rows):\nG551D\nF508del\np.Gly551Asp\nE1371Q\nM470V`}
              value={textInput}
              onChange={(e) => setTextInput(e.target.value)}
            />
            <div style={{ display: 'flex', justifyContent: 'space-between', marginTop: '6px', fontSize: '12px', color: '#64748b' }}>
              <span>
                {currentIdentifiers.length} / 500 non-empty variants entered
              </span>
              <span>Supported: Canonical protein changes (G551D), HGVS (p.Gly551Asp), legacy names</span>
            </div>
          </div>
        ) : (
          <div
            style={{
              padding: '24px',
              border: '2px dashed #cbd5e1',
              borderRadius: '8px',
              textAlign: 'center',
              background: '#f8fafc',
            }}
          >
            <Upload style={{ width: '32px', height: '32px', color: '#0f766e', margin: '0 auto 8px' }} />
            <p style={{ margin: '0 0 6px', fontWeight: 600 }}>Select a CSV file containing variant identifiers</p>
            <p style={{ margin: '0 0 12px', fontSize: '12px', color: '#64748b' }}>
              Limit 2 MB (up to 500 rows). Column header auto-detects `variant`, `identifier`, or uses the 1st column.
            </p>
            <input
              type="file"
              accept=".csv"
              onChange={handleFileUpload}
              style={{ display: 'inline-block', fontSize: '13px' }}
            />
            {csvFileName && (
              <p style={{ marginTop: '8px', fontSize: '12px', color: '#0f766e', fontWeight: 600 }}>
                Loaded: {csvFileName} ({csvIdentifiers.length} identifiers extracted)
              </p>
            )}
            {csvError && (
              <p style={{ marginTop: '8px', fontSize: '12px', color: '#b91c1c' }}>
                {csvError}
              </p>
            )}
          </div>
        )}

        {error && (
          <div
            style={{
              marginTop: '12px',
              padding: '8px 12px',
              background: '#fef2f2',
              color: '#991b1b',
              borderRadius: '6px',
              fontSize: '13px',
              display: 'flex',
              alignItems: 'center',
              gap: '8px',
            }}
          >
            <AlertTriangle style={{ width: '16px' }} />
            <span>{error}</span>
          </div>
        )}

        <div style={{ marginTop: '16px', display: 'flex', justifyContent: 'flex-end', gap: '8px' }}>
          <button
            className="primary"
            style={{ padding: '8px 20px', fontSize: '14px', fontWeight: 600 }}
            onClick={handleRunBatch}
            disabled={loading || currentIdentifiers.length === 0}
          >
            {loading ? 'Analyzing Batch...' : `Run Batch Analysis (${currentIdentifiers.length})`}
          </button>
        </div>
      </article>

      {/* Summary Metrics */}
      {response && (
        <>
          <div className="metric-strip" style={{ marginBottom: '16px' }}>
            <article className="metric">
              <div className="metric-label">Input Rows</div>
              <strong>{response.total_input_rows.toLocaleString()}</strong>
              <p>Submitted lines</p>
            </article>
            <article className="metric">
              <div className="metric-label">Unique Keys</div>
              <strong>{response.unique_identifiers.toLocaleString()}</strong>
              <p>Deduplicated queries</p>
            </article>
            <article className="metric">
              <div className="metric-label">Resolved</div>
              <strong style={{ color: '#0f766e' }}>{response.resolved_count.toLocaleString()}</strong>
              <p>Found in database</p>
            </article>
            <article className="metric">
              <div className="metric-label">Unresolved</div>
              <strong style={{ color: response.unresolved_count > 0 ? '#b91c1c' : '#64748b' }}>
                {response.unresolved_count.toLocaleString()}
              </strong>
              <p>Missing or syntax error</p>
            </article>
            <article className="metric">
              <div className="metric-label">Structural Evidence</div>
              <strong style={{ color: '#0f766e' }}>{response.evidence_available_count.toLocaleString()}</strong>
              <p>Mapped to Cryo-EM models</p>
            </article>
          </div>

          {/* Results Table Section */}
          <article className="card">
            <div
              style={{
                display: 'flex',
                justifyContent: 'space-between',
                alignItems: 'center',
                flexWrap: 'wrap',
                gap: '12px',
                marginBottom: '12px',
              }}
            >
              <div>
                <h3 style={{ margin: 0, fontSize: '16px' }}>Batch Results</h3>
                <small style={{ color: '#64748b' }}>
                  Batch ID: {response.batch_id} · Analyzed at {new Date(response.timestamp).toLocaleTimeString()}
                </small>
              </div>

              <div style={{ display: 'flex', gap: '8px', alignItems: 'center' }}>
                <div style={{ position: 'relative' }}>
                  <input
                    type="text"
                    placeholder="Search results..."
                    value={tableSearch}
                    onChange={(e) => {
                      setTableSearch(e.target.value);
                      setPage(0);
                    }}
                    style={{
                      padding: '5px 8px 5px 28px',
                      fontSize: '12px',
                      borderRadius: '4px',
                      border: '1px solid #cbd5e1',
                    }}
                  />
                  <Search
                    style={{
                      position: 'absolute',
                      left: '8px',
                      top: '7px',
                      width: '13px',
                      color: '#94a3b8',
                    }}
                  />
                </div>

                <select
                  value={statusFilter}
                  onChange={(e) => {
                    setStatusFilter(e.target.value as any);
                    setPage(0);
                  }}
                  style={{
                    padding: '5px 8px',
                    fontSize: '12px',
                    borderRadius: '4px',
                    border: '1px solid #cbd5e1',
                  }}
                >
                  <option value="ALL">All Statuses ({response.results.length})</option>
                  <option value="RESOLVED">Resolved Only ({response.resolved_count})</option>
                  <option value="UNRESOLVED">Unresolved Only ({response.unresolved_count})</option>
                </select>

                <button
                  className="secondary"
                  style={{ padding: '6px 12px', fontSize: '12px' }}
                  onClick={handleExportResults}
                >
                  <Download style={{ width: '13px', marginRight: '4px', display: 'inline-block' }} />
                  Export CSV
                </button>
              </div>
            </div>

            <div className="table-wrap">
              <table>
                <thead>
                  <tr>
                    <th>#</th>
                    <th>Input identifier</th>
                    <th>Normalized variant</th>
                    <th>Status</th>
                    <th>Sources & Determinations</th>
                    <th>Structure coverage</th>
                    <th>Warnings / Caveats</th>
                    <th>Actions</th>
                  </tr>
                </thead>
                <tbody>
                  {paginatedRows.map((r) => {
                    const s5 = r.structures?.['5UAK'];
                    const s6 = r.structures?.['6MSM'];
                    const s8 = r.structures?.['8EIQ'];

                    return (
                      <tr key={r.row_number}>
                        <td style={{ color: '#64748b' }}>{r.row_number}</td>
                        <td>
                          <b>{r.original_input || <i style={{ color: '#94a3b8' }}>[blank]</i>}</b>
                        </td>
                        <td>
                          {r.normalized_identifier ? (
                            <span>
                              <b>{r.normalized_identifier}</b>
                              {r.canonical_position && (
                                <small style={{ display: 'block', color: '#64748b' }}>
                                  Pos {r.canonical_position} · {r.domain_name ?? 'Domain —'}
                                </small>
                              )}
                            </span>
                          ) : (
                            <span style={{ color: '#94a3b8' }}>—</span>
                          )}
                        </td>
                        <td>
                          <span
                            className={`badge ${
                              r.status === 'RESOLVED'
                                ? 'good'
                                : r.status === 'NOT_FOUND'
                                  ? 'warn'
                                  : 'neutral'
                            }`}
                          >
                            {r.status}
                          </span>
                        </td>
                        <td>
                          {r.source_record_count > 0 ? (
                            <div>
                              <span>{r.source_record_count} record(s)</span>
                              {r.source_determinations.length > 0 && (
                                <small style={{ display: 'block', color: '#475569' }}>
                                  {r.source_determinations.join('; ')}
                                </small>
                              )}
                            </div>
                          ) : (
                            <span style={{ color: '#94a3b8' }}>None</span>
                          )}
                        </td>
                        <td>
                          {r.has_structural_evidence ? (
                            <div style={{ display: 'flex', gap: '4px' }}>
                              <span
                                className={`badge ${
                                  s5?.mapping_status === 'MAPPED_OBSERVED' ? 'teal' : 'neutral'
                                }`}
                                title={`5UAK: ${s5?.mapping_status}`}
                              >
                                5UAK
                              </span>
                              <span
                                className={`badge ${
                                  s6?.mapping_status === 'MAPPED_OBSERVED' ? 'teal' : 'neutral'
                                }`}
                                title={`6MSM: ${s6?.mapping_status}`}
                              >
                                6MSM
                              </span>
                              <span
                                className={`badge ${
                                  s8?.mapping_status === 'MAPPED_OBSERVED'
                                    ? 'teal'
                                    : s8?.failure_reason === 'CONSTRUCT_DELETION'
                                      ? 'warn'
                                      : 'neutral'
                                }`}
                                title={`8EIQ: ${s8?.mapping_status ?? s8?.failure_reason}`}
                              >
                                8EIQ
                              </span>
                            </div>
                          ) : (
                            <span style={{ color: '#94a3b8' }}>—</span>
                          )}
                        </td>
                        <td>
                          {r.status !== 'RESOLVED' ? (
                            <small style={{ color: '#b91c1c' }}>{r.failure_reason}</small>
                          ) : r.warnings.length > 0 ? (
                            <span
                              className="badge warn"
                              title={r.warnings.map((w) => w.message).join(' | ')}
                            >
                              {r.warnings.length} warning(s)
                            </span>
                          ) : (
                            <span style={{ color: '#0f766e', fontSize: '11px' }}>✓ Normal</span>
                          )}
                        </td>
                        <td>
                          {r.status === 'RESOLVED' && (
                            <div style={{ display: 'flex', gap: '4px' }}>
                              <button
                                className="icon-button"
                                title="Inspect in 3D Structural Laboratory"
                                onClick={() =>
                                  onInspect3D(
                                    r.normalized_identifier ?? r.original_input,
                                    r.canonical_position ?? undefined
                                  )
                                }
                              >
                                <Activity style={{ width: '13px' }} />
                              </button>
                              <button
                                className="icon-button"
                                title="Open Structural Evidence Detail"
                                onClick={() =>
                                  onSelectVariant(
                                    r.normalized_identifier ?? r.original_input
                                  )
                                }
                              >
                                <ArrowRight style={{ width: '13px' }} />
                              </button>
                            </div>
                          )}
                        </td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </div>

            {/* Pagination */}
            {totalPages > 1 && (
              <div className="pagination" style={{ marginTop: '12px' }}>
                <span>
                  Showing {page * pageSize + 1}–
                  {Math.min((page + 1) * pageSize, filteredResults.length)} of{' '}
                  {filteredResults.length}
                </span>
                <div>
                  <button
                    className="icon-button"
                    disabled={page === 0}
                    onClick={() => setPage((p) => p - 1)}
                    aria-label="Previous page"
                  >
                    <ChevronLeft />
                  </button>
                  <button
                    className="icon-button"
                    disabled={page >= totalPages - 1}
                    onClick={() => setPage((p) => p + 1)}
                    aria-label="Next page"
                  >
                    <ChevronRight />
                  </button>
                </div>
              </div>
            )}
          </article>
        </>
      )}
    </>
  );
}
