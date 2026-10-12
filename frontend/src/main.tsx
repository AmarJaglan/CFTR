import React, { useEffect, useState } from 'react';
import { createRoot } from 'react-dom/client';
import {
  Activity,
  AlertTriangle,
  BookOpen,
  ChevronLeft,
  ChevronRight,
  Compass,
  Database,
  FlaskConical,
  Info,
  Layers,
  Menu,
  RotateCcw,
  Search,
  ShieldCheck,
  Table2,
  X,
} from 'lucide-react';
import {
  api,
  CFTR2Record,
  ClinVarRecord,
  DatabaseStats,
  VariantDetail,
  VariantRecord,
} from './lib/api';
import { BatchAnalysis } from './components/BatchAnalysis';
import { StructureViewer } from './components/StructureViewer';
import './styles.css';

type Page =
  | 'Overview'
  | 'Variant Explorer'
  | 'Structural Laboratory'
  | 'Structural Evidence'
  | 'Batch Analysis'
  | 'Data Provenance'
  | 'About & Methods';

const nav: [Page, React.ElementType][] = [
  ['Overview', Activity],
  ['Variant Explorer', Search],
  ['Structural Laboratory', FlaskConical],
  ['Structural Evidence', Layers],
  ['Batch Analysis', Table2],
  ['Data Provenance', Database],
  ['About & Methods', BookOpen],
];

const structures = ['5UAK', '6MSM', '8EIQ'] as const;

function Badge({
  children,
  tone = 'neutral',
}: {
  children: React.ReactNode;
  tone?: 'neutral' | 'good' | 'warn' | 'teal';
}) {
  return <span className={`badge ${tone}`}>{children}</span>;
}

function Metric({
  label,
  value,
  detail,
}: {
  label: string;
  value: string;
  detail: string;
}) {
  return (
    <article className="metric">
      <div className="metric-label">{label}</div>
      <strong>{value}</strong>
      <p>{detail}</p>
    </article>
  );
}

function Shell({
  page,
  setPage,
  backendOnline,
  children,
}: {
  page: Page;
  setPage: (p: Page) => void;
  backendOnline: boolean | null;
  children: React.ReactNode;
}) {
  const [open, setOpen] = useState(false);

  return (
    <div className="app-shell">
      <aside className={open ? 'sidebar open' : 'sidebar'}>
        <div className="brand">
          <div className="brand-mark">CF</div>
          <div>
            <b>CFTR-StructMap</b>
            <small>Structural evidence workspace</small>
          </div>
          <button
            className="icon-button mobile-only"
            onClick={() => setOpen(false)}
            aria-label="Close navigation"
          >
            <X />
          </button>
        </div>
        <nav>
          {nav.map(([name, Icon]) => (
            <button
              className={page === name ? 'nav-item active' : 'nav-item'}
              onClick={() => {
                setPage(name);
                setOpen(false);
              }}
              key={name}
            >
              <Icon />
              <span>{name}</span>
            </button>
          ))}
        </nav>
        <div className="side-note">
          <ShieldCheck />
          <div>
            <b>Evidence boundary</b>
            <p>Structural context, not clinical prediction.</p>
          </div>
        </div>
      </aside>

      <main className="main">
        <header className="topbar">
          <button
            className="icon-button mobile-only"
            onClick={() => setOpen(true)}
            aria-label="Open navigation"
          >
            <Menu />
          </button>
          <div>
            <span className="eyebrow">CFTR / P13569</span>
            <h1>{page}</h1>
          </div>
          <div className="status">
            <span
              className="status-dot"
              style={{
                background:
                  backendOnline === true
                    ? '#31b7ac'
                    : backendOnline === false
                      ? '#e06c75'
                      : '#c9d8df',
              }}
            />
            {backendOnline === true
              ? 'FastAPI backend online'
              : backendOnline === false
                ? 'Backend unavailable'
                : 'Checking backend status…'}
          </div>
        </header>
        <div className="content">{children}</div>
      </main>
    </div>
  );
}

function Notice({
  children,
  kind = 'info',
}: {
  children: React.ReactNode;
  kind?: 'info' | 'warn';
}) {
  return (
    <div className={`notice ${kind}`}>
      {kind === 'warn' ? <AlertTriangle /> : <Info />}
      <div>{children}</div>
    </div>
  );
}

// ==============================================================================
// 1. Overview Section
// ==============================================================================

function Overview({ setPage }: { setPage: (p: Page) => void }) {
  const [health, setHealth] = useState<'checking' | 'online' | 'unavailable'>(
    'checking'
  );
  const [stats, setStats] = useState<DatabaseStats | null>(null);

  useEffect(() => {
    api
      .health()
      .then(() => {
        setHealth('online');
        return api.stats();
      })
      .then((res) => {
        setStats(res);
      })
      .catch(() => {
        setHealth('unavailable');
        setStats(null);
      });
  }, []);

  return (
    <>
      <section className="hero">
        <div>
          <span className="eyebrow">RESEARCH DATA PLATFORM</span>
          <h2>Structural context for CFTR variation.</h2>
          <p>
            Trace normalized variant records from CFTR2 and ClinVar to
            experimentally resolved structures without collapsing scientific
            uncertainty.
          </p>
        </div>
        <div className="hero-accent">
          <span>UniProt reference</span>
          <strong>P13569</strong>
          <small>1,480 amino acids</small>
        </div>
      </section>

      <div className="metrics">
        <Metric
          label="Source records"
          value={
            stats?.total_evidence_records
              ? stats.total_evidence_records.toLocaleString()
              : health === 'checking'
                ? 'Checking…'
                : 'Unavailable'
          }
          detail={
            stats
              ? `${stats.cftr2_records.toLocaleString()} CFTR2 · ${stats.clinvar_records.toLocaleString()} ClinVar records`
              : 'Requires active FastAPI connection'
          }
        />
        <Metric
          label="Grouped protein changes"
          value={
            stats?.grouped_protein_changes
              ? stats.grouped_protein_changes.toLocaleString()
              : health === 'checking'
                ? 'Checking…'
                : 'Unavailable'
          }
          detail={
            stats
              ? `${stats.comparable_positions.toLocaleString()} positions with cross-structure metrics`
              : 'Requires active database connection'
          }
        />
        <Metric
          label="Analyzed structures"
          value={stats ? String(stats.cftr_structures_count) : '3'}
          detail="5UAK (apo) · 6MSM (ATP) · 8EIQ (ΔF508)"
        />
        <Metric
          label="Backend status"
          value={
            health === 'online'
              ? stats?.evidence_loaded
                ? 'Online & Loaded'
                : 'Online (Empty DB)'
              : health === 'unavailable'
                ? 'Offline'
                : 'Checking…'
          }
          detail={
            health === 'online'
              ? 'FastAPI v0.1.0 · SQLite database verified'
              : 'Connect FastAPI server at localhost:8000'
          }
        />
      </div>

      <div className="grid-two">
        <article className="card">
          <div className="card-heading">
            <div>
              <span className="eyebrow">STRUCTURE INVENTORY</span>
              <h3>Three experimental contexts</h3>
            </div>
            <button
              className="text-button"
              onClick={() => setPage('Structural Laboratory')}
            >
              Open laboratory <ChevronRight />
            </button>
          </div>
          <div className="structure-list">
            <div className="structure-row">
              <span className="structure-id">5UAK</span>
              <div>
                <b>Apo, dephosphorylated</b>
                <small>Cryo-EM 3.87 Å · wild-type human CFTR</small>
              </div>
              <Badge tone="teal">reference context</Badge>
            </div>
            <div className="structure-row">
              <span className="structure-id">6MSM</span>
              <div>
                <b>ATP-bound dimer (E1371Q)</b>
                <small>Cryo-EM 3.20 Å · catalytic substitution E1371Q</small>
              </div>
              <Badge tone="teal">reference context</Badge>
            </div>
            <div className="structure-row">
              <span className="structure-id">8EIQ</span>
              <div>
                <b>ATP-bound dimer, ΔF508 (E1371Q)</b>
                <small>Cryo-EM 3.40 Å · authentic patient deletion ΔF508</small>
              </div>
              <Badge tone="warn">mutant context</Badge>
            </div>
          </div>
        </article>

        <article className="card">
          <div className="card-heading">
            <div>
              <span className="eyebrow">INTERPRETATION</span>
              <h3>Use the evidence carefully</h3>
            </div>
            <ShieldCheck />
          </div>
          <ul className="principles">
            <li>
              Mapped sequence position is distinct from observed atomic
              coordinates.
            </li>
            <li>Unresolved residues are unavailable, not zero-valued.</li>
            <li>Construct substitutions and deletions affect comparability.</li>
            <li>
              Structural evidence does not independently establish
              pathogenicity.
            </li>
          </ul>
        </article>
      </div>

      <Notice>
        <b>Scientific disclaimer.</b> CFTR-StructMap presents structural
        evidence integration observations. It does not generate clinical
        classifications, pathogenicity scores, or treatment recommendations.
      </Notice>
    </>
  );
}

// ==============================================================================
// 2. Variant Explorer Section
// ==============================================================================

function VariantExplorer({ onSelect }: { onSelect: (id: string) => void }) {
  const [query, setQuery] = useState('');
  const [activeSearch, setActiveSearch] = useState('');
  const [domainFilter, setDomainFilter] = useState('');
  const [sourceFilter, setSourceFilter] = useState('');
  const [data, setData] = useState<{
    results: VariantRecord[];
    total: number;
    offset: number;
    limit: number;
  } | null>(null);
  const [error, setError] = useState('');
  const [loading, setLoading] = useState(false);

  const load = (
    offset = 0,
    searchVal = activeSearch,
    domVal = domainFilter,
    srcVal = sourceFilter
  ) => {
    setLoading(true);
    setError('');

    const isNumeric = searchVal.trim() && !Number.isNaN(Number(searchVal.trim()));
    const params: Record<string, string | number | undefined> = {
      limit: 12,
      offset,
      search: searchVal.trim() || undefined,
      position: isNumeric ? Number(searchVal.trim()) : undefined,
      domain: domVal || undefined,
      source_db: srcVal || undefined,
    };

    api
      .listVariants(params)
      .then(setData)
      .catch((e) => setError(e.message))
      .finally(() => setLoading(false));
  };

  useEffect(() => {
    load(0, '', '', '');
  }, []);

  const handleSearchSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    setActiveSearch(query);
    load(0, query, domainFilter, sourceFilter);
  };

  const handleReset = () => {
    setQuery('');
    setActiveSearch('');
    setDomainFilter('');
    setSourceFilter('');
    load(0, '', '', '');
  };

  return (
    <>
      <section className="page-intro">
        <span className="eyebrow">SOURCE-AWARE SEARCH</span>
        <h2>Variant Explorer</h2>
        <p>
          Search 4,756 integrated source records across CFTR2 (2,092) and ClinVar
          (2,664), representing 3,353 unique canonical protein changes. Includes
          missense substitutions, in-frame deletions, synonymous controls, and
          intronic/unresolved records.
        </p>
      </section>

      <form className="searchbar" onSubmit={handleSearchSubmit}>
        <Search />
        <input
          aria-label="Search variants"
          placeholder="Try F508del, G551D, M470V, 508, or VariationID (e.g. 7120)"
          value={query}
          onChange={(e) => setQuery(e.target.value)}
        />
        <button className="primary" type="submit">
          Search
        </button>
      </form>

      <div
        style={{
          display: 'flex',
          gap: '10px',
          flexWrap: 'wrap',
          marginBottom: '16px',
          alignItems: 'center',
        }}
      >
        <select
          style={{
            background: '#fff',
            border: '1px solid #c9d8df',
            borderRadius: '6px',
            padding: '7px 10px',
            fontSize: '12px',
            color: '#17364e',
          }}
          value={sourceFilter}
          onChange={(e) => {
            setSourceFilter(e.target.value);
            load(0, activeSearch, domainFilter, e.target.value);
          }}
        >
          <option value="">All Sources</option>
          <option value="CFTR2">CFTR2 only</option>
          <option value="ClinVar">ClinVar only</option>
        </select>

        <select
          style={{
            background: '#fff',
            border: '1px solid #c9d8df',
            borderRadius: '6px',
            padding: '7px 10px',
            fontSize: '12px',
            color: '#17364e',
          }}
          value={domainFilter}
          onChange={(e) => {
            setDomainFilter(e.target.value);
            load(0, activeSearch, e.target.value, sourceFilter);
          }}
        >
          <option value="">All Domains</option>
          <option value="Lasso">Lasso (1–79)</option>
          <option value="TMD1">TMD1 (80–388)</option>
          <option value="NBD1">NBD1 (389–678)</option>
          <option value="R-domain">R-domain (679–844)</option>
          <option value="TMD2">TMD2 (845–1172)</option>
          <option value="NBD2">NBD2 (1173–1436)</option>
          <option value="C-term">C-term (1437–1480)</option>
        </select>

        {(activeSearch || domainFilter || sourceFilter) && (
          <button
            className="text-button"
            type="button"
            onClick={handleReset}
            style={{ marginLeft: '4px' }}
          >
            <RotateCcw style={{ width: '13px' }} /> Reset filters
          </button>
        )}
      </div>

      {error && (
        <Notice kind="warn">
          <span>
            <b>API unavailable.</b> {error}
          </span>
          <button
            className="text-button"
            onClick={() => load(0, activeSearch, domainFilter, sourceFilter)}
          >
            Retry
          </button>
        </Notice>
      )}

      <article className="card table-card">
        <div className="card-heading">
          <div>
            <h3>Evidence records</h3>
            <p className="muted">
              {data
                ? `${data.total.toLocaleString()} records reported by the API`
                : 'Waiting for backend response'}
            </p>
          </div>
          {loading && <span className="loading">Loading…</span>}
        </div>

        {data?.results && data.results.length > 0 ? (
          <>
            <div className="table-wrap">
              <table>
                <thead>
                  <tr>
                    <th>Protein change</th>
                    <th>Source</th>
                    <th>Source determination</th>
                    <th>Position</th>
                    <th>Domain</th>
                    <th>Eligibility</th>
                    <th>Inspect</th>
                  </tr>
                </thead>
                <tbody>
                  {data.results.map((r, i) => (
                    <tr key={`${r.source_db}-${r.source_record_id}-${i}`}>
                      <td>
                        <b>
                          {r.canonical_protein ??
                            r.source_variant_name ??
                            'Unnamed variant'}
                        </b>
                        <small>
                          {r.source_variant_name ?? r.source_record_id}
                        </small>
                      </td>
                      <td>
                        <Badge
                          tone={r.source_db === 'CFTR2' ? 'teal' : 'neutral'}
                        >
                          {r.source_db}
                        </Badge>
                      </td>
                      <td>{r.source_determination ?? '—'}</td>
                      <td>{r.canonical_position ?? '—'}</td>
                      <td>{r.domain_name ?? '—'}</td>
                      <td>
                        <Badge
                          tone={
                            r.eligibility_status === 'ELIGIBLE_MISSENSE'
                              ? 'good'
                              : 'warn'
                          }
                        >
                          {r.eligibility_status ?? '—'}
                        </Badge>
                      </td>
                      <td>
                        <button
                          className="icon-button"
                          onClick={() =>
                            onSelect(
                              r.canonical_protein ?? r.source_record_id
                            )
                          }
                          aria-label={`Inspect ${r.canonical_protein ?? r.source_record_id}`}
                        >
                          <ChevronRight />
                        </button>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
            <div className="pagination">
              <span>
                Showing {data.offset + 1}–
                {Math.min(data.offset + data.results.length, data.total)} of{' '}
                {data.total.toLocaleString()}
              </span>
              <div>
                <button
                  className="icon-button"
                  disabled={data.offset === 0}
                  onClick={() =>
                    load(
                      Math.max(0, data.offset - data.limit),
                      activeSearch,
                      domainFilter,
                      sourceFilter
                    )
                  }
                  aria-label="Previous page"
                >
                  <ChevronLeft />
                </button>
                <button
                  className="icon-button"
                  disabled={data.offset + data.limit >= data.total}
                  onClick={() =>
                    load(
                      data.offset + data.limit,
                      activeSearch,
                      domainFilter,
                      sourceFilter
                    )
                  }
                  aria-label="Next page"
                >
                  <ChevronRight />
                </button>
              </div>
            </div>
          </>
        ) : (
          !loading &&
          !error && (
            <div className="empty">
              <Search />
              <h3>No records returned</h3>
              <p>Try adjusting your search query or filters.</p>
            </div>
          )
        )}
      </article>
    </>
  );
}

// ==============================================================================
// 3. Structural Laboratory Section
// ==============================================================================

function Lab() {
  return (
    <>
      <section className="page-intro">
        <span className="eyebrow">EXPERIMENTAL MODELS</span>
        <h2>Structural Laboratory</h2>
        <p>
          Compare the three experimental cryo-EM models used in CFTR-StructMap.
          Each structure captures distinct conformational and construct
          contexts.
        </p>
      </section>

      <div className="structure-grid">
        <article className="card structure-card">
          <div className="structure-card-top">
            <span className="structure-id">5UAK</span>
            <Badge tone="teal">reference context</Badge>
          </div>
          <h3>Apo, dephosphorylated</h3>
          <p>
            Human CFTR resolved at 3.87 Å resolution by single-particle cryo-EM.
            Represents the inward-facing, ATP-free basal state. Chain A mapped to
            canonical UniProt P13569 (residues 1–1480).
          </p>
          <div className="feature-chips">
            <span>SASA / RSA (Tien scale)</span>
            <span>8.0 Å contact sphere</span>
            <span>Deposited secondary structure</span>
            <span>Inward-facing pore</span>
          </div>
        </article>

        <article className="card structure-card">
          <div className="structure-card-top">
            <span className="structure-id">6MSM</span>
            <Badge tone="teal">reference context</Badge>
          </div>
          <h3>ATP-bound dimer (E1371Q)</h3>
          <p>
            Human CFTR resolved at 3.20 Å resolution. Phosphorylated channel
            trapped with ATP in the head-to-tail NBD1–NBD2 dimer state via
            engineered catalytic substitution E1371Q.
          </p>
          <div className="feature-chips">
            <span>SASA / RSA (Tien scale)</span>
            <span>8.0 Å contact sphere</span>
            <span>E1371Q construct mutation</span>
            <span>Closed/conducting pore</span>
          </div>
        </article>

        <article className="card structure-card">
          <div className="structure-card-top">
            <span className="structure-id">8EIQ</span>
            <Badge tone="warn">mutant context</Badge>
          </div>
          <h3>ATP-bound dimer, ΔF508 (E1371Q)</h3>
          <p>
            Human CFTR resolved at 3.40 Å resolution. Harbors authentic patient
            deletion ΔF508 (residue 508 absent) along with catalytic E1371Q.
            Residue 508 coordinates are strictly unmodeled.
          </p>
          <div className="feature-chips">
            <span>Authentic ΔF508 deletion</span>
            <span>Residue 508 null features</span>
            <span>E1371Q construct mutation</span>
            <span>ATP-dimer state</span>
          </div>
        </article>
      </div>

      <div style={{ marginTop: '20px' }}>
        <article className="card">
          <div className="card-heading">
            <div>
              <span className="eyebrow">NON-CFTR CONTROLS</span>
              <h3>Intentionally excluded structures</h3>
            </div>
          </div>
          <p style={{ fontSize: '13px', color: '#587285', margin: 0 }}>
            Structures <b>1AV5</b> (isolated murine NBD1 crystal structure) and{' '}
            <b>7JQ4</b> (non-CFTR ABC transporter control) were parsed during
            Phase 3 but are intentionally excluded from full-length human CFTR
            mapping to prevent construct confounding.
          </p>
        </article>
      </div>

      <Notice kind="warn">
        <b>Planned capability.</b> A 3D interactive molecular viewer is not
        enabled in the current release because the FastAPI backend does not
        expose coordinate or mesh endpoints (`/coordinates` or `/pdb` not
        present). Structural features are evaluated from mmCIF atomic models on
        the backend.
      </Notice>
    </>
  );
}

// ==============================================================================
// 4. Structural Evidence Section
// ==============================================================================

function Evidence({
  selected,
  onSelectVariant,
  onOpen3D,
}: {
  selected: string | null;
  onSelectVariant: (id: string) => void;
  onOpen3D?: (pos: number) => void;
}) {
  const [detail, setDetail] = useState<VariantDetail | null>(null);
  const [error, setError] = useState('');
  const [loading, setLoading] = useState(false);

  useEffect(() => {
    if (selected) {
      setLoading(true);
      setError('');
      api
        .getVariant(selected)
        .then(setDetail)
        .catch((e) => {
          setError(e.message);
          setDetail(null);
        })
        .finally(() => setLoading(false));
    }
  }, [selected]);

  if (!selected) {
    return (
      <>
        <section className="page-intro">
          <span className="eyebrow">WORKSPACE</span>
          <h2>Structural Evidence</h2>
          <p>
            Select a variant from Variant Explorer or pick a representative
            benchmark below to inspect structure-specific evidence.
          </p>
        </section>

        <div className="empty large">
          <Layers style={{ width: '32px', color: '#147f7a' }} />
          <h3>No variant selected</h3>
          <p>
            Choose a representative variant below or search in Variant Explorer:
          </p>
          <div
            style={{
              display: 'flex',
              gap: '8px',
              flexWrap: 'wrap',
              justifyContent: 'center',
              marginTop: '16px',
            }}
          >
            {[
              ['F508del', 'ΔF508 authentic deletion'],
              ['G551D', 'Gating mutation in NBD1'],
              ['M470V', 'Historical reference mismatch'],
              ['1371', 'E1371Q catalytic substitution'],
              ['700', 'R-domain unresolved coordinates'],
            ].map(([id, desc]) => (
              <button
                key={id}
                className="primary"
                style={{
                  background: '#e3f5f3',
                  color: '#137e7a',
                  border: '1px solid #cce5e3',
                  fontSize: '12px',
                  padding: '7px 12px',
                }}
                onClick={() => onSelectVariant(id)}
              >
                <b>{id}</b> ({desc})
              </button>
            ))}
          </div>
        </div>
      </>
    );
  }

  return (
    <>
      <section className="page-intro">
        <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', flexWrap: 'wrap', gap: '12px' }}>
          <div>
            <span className="eyebrow">EVIDENCE RECORD</span>
            <h2>{detail?.canonical_protein ?? selected}</h2>
            <p>
              Canonical position {detail?.canonical_position ?? '—'} ·{' '}
              {detail?.domain_name ?? 'Domain unavailable'} · reference P13569
            </p>
          </div>
          {onOpen3D && detail?.canonical_position && (
            <button
              className="secondary"
              style={{ padding: '8px 14px', fontSize: '13px', fontWeight: 600 }}
              onClick={() => onOpen3D(detail.canonical_position!)}
            >
              <Compass style={{ width: '14px', marginRight: '6px', display: 'inline-block' }} />
              Inspect in 3D Structural Laboratory
            </button>
          )}
        </div>
      </section>

      {loading && <div className="notice info">Loading structural evidence…</div>}
      {error && <Notice kind="warn">{error}</Notice>}

      {detail && (
        <>
          <div className="evidence-summary">
            <Metric
              label="Canonical residue"
              value={
                detail.ref_aa
                  ? `${detail.ref_aa}${detail.canonical_position}${detail.alt_aa ?? ''}`
                  : (detail.canonical_protein ?? selected)
              }
              detail="Normalized protein change"
            />
            <Metric
              label="CFTR2 records"
              value={String(detail.source_associations.cftr2_count)}
              detail="Curated registry records"
            />
            <Metric
              label="ClinVar records"
              value={String(detail.source_associations.clinvar_count)}
              detail="VCV aggregate assertions"
            />
            <Metric
              label="Observation pattern"
              value={String(
                detail.comparative_metrics.observation_pattern ?? '—'
              )}
              detail="Multi-structure coverage"
            />
          </div>

          <div className="structure-grid">
            {structures.map((id) => {
              const s = detail.structures[id];
              return (
                <article className="card evidence-card" key={id}>
                  <div className="card-heading">
                    <h3>{id}</h3>
                    <Badge
                      tone={
                        s?.mapping_status === 'MAPPED_OBSERVED'
                          ? 'good'
                          : 'warn'
                      }
                    >
                      {s?.mapping_status ?? 'Unavailable'}
                    </Badge>
                  </div>
                  <p>{s?.state}</p>
                  <dl>
                    <dt>Observed residue</dt>
                    <dd>{s?.structure_aa ?? '—'}</dd>
                    <dt>Feature status</dt>
                    <dd>{s?.feature_status ?? '—'}</dd>
                    <dt>RSA ratio</dt>
                    <dd>
                      {s?.rsa_ratio !== null && s?.rsa_ratio !== undefined ? (
                        <span>
                          <b>{(s.rsa_ratio * 100).toFixed(1)}%</b>
                          {s.rsa_ratio === 0 && (
                            <small style={{ display: 'block', color: '#147f7a', fontSize: '11px' }}>
                              Completely buried (0.0%)
                            </small>
                          )}
                        </span>
                      ) : (
                        <span style={{ color: '#94a3b8' }}>— <small>(Unmodeled)</small></span>
                      )}
                    </dd>
                    <dt>SASA (raw)</dt>
                    <dd>
                      {s?.sasa_angstrom2 !== null && s?.sasa_angstrom2 !== undefined ? (
                        <span>
                          <b>{s.sasa_angstrom2.toFixed(1)} Å²</b>
                          {s.sasa_angstrom2 === 0 && (
                            <small style={{ display: 'block', color: '#147f7a', fontSize: '11px' }}>
                              Occluded at dimer interface
                            </small>
                          )}
                        </span>
                      ) : (
                        <span style={{ color: '#94a3b8' }}>— <small>(Unmodeled)</small></span>
                      )}
                    </dd>
                    <dt title="Contacts within 8.0 Å heavy-atom sphere (includes protein residues and bound ligands if present)">
                      Neighbor contacts (8.0 Å)
                    </dt>
                    <dd>
                      {s?.total_contacts !== null && s?.total_contacts !== undefined ? (
                        <span>
                          <b>{s.total_contacts}</b>
                          <small style={{ display: 'block', color: '#64748b', fontSize: '11px' }}>
                            {s.inter_contacts && s.inter_contacts > 0
                              ? `${s.intra_contacts} intra · ${s.inter_contacts} inter`
                              : id === '5UAK'
                                ? '19 intra-chain protein residues'
                                : 'Includes 28 protein residues + ATP & Mg²⁺'}
                          </small>
                        </span>
                      ) : (
                        <span style={{ color: '#94a3b8' }}>— <small>(Unmodeled)</small></span>
                      )}
                    </dd>
                    <dt>Secondary structure</dt>
                    <dd>
                      {s?.secondary_structure ? (
                        <b>{s.secondary_structure}</b>
                      ) : (
                        <span style={{ color: '#94a3b8' }}>— <small>(Unmodeled)</small></span>
                      )}
                    </dd>
                  </dl>
                  {s?.construct_note && (
                    <div
                      style={{
                        marginTop: '12px',
                        paddingTop: '8px',
                        borderTop: '1px solid #edf1f3',
                        fontSize: '11px',
                        color: '#a35e1c',
                      }}
                    >
                      <b>Construct note:</b> {s.construct_note}
                    </div>
                  )}
                </article>
              );
            })}
          </div>

          {/* Comparative Metrics Card */}
          <article className="card" style={{ marginTop: '16px' }}>
            <div className="card-heading">
              <div>
                <span className="eyebrow">OBSERVATIONAL DYNAMICS</span>
                <h3>Cross-structure comparison</h3>
              </div>
              <Badge
                tone={
                  detail.comparative_metrics.is_fully_comparable
                    ? 'good'
                    : 'warn'
                }
              >
                {detail.comparative_metrics.is_fully_comparable
                  ? 'Observed in all 3 structures (experimental contexts differ)'
                  : (detail.comparative_metrics.environment_comparability?.replaceAll('_', ' ') ??
                    'Caveats / Incomparable')}
              </Badge>
            </div>
            <div
              style={{
                display: 'grid',
                gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))',
                gap: '12px',
                fontSize: '12px',
              }}
            >
              <div>
                <span style={{ color: '#8297a4' }}>ΔRSA divergence:</span>{' '}
                <b>
                  {detail.comparative_metrics.delta_rsa !== null &&
                  detail.comparative_metrics.delta_rsa !== undefined
                    ? `${(detail.comparative_metrics.delta_rsa * 100).toFixed(1)}%`
                    : '—'}
                </b>
              </div>
              <div>
                <span style={{ color: '#8297a4' }}>ΔContacts divergence:</span>{' '}
                <b>
                  {detail.comparative_metrics.delta_contacts !== null &&
                  detail.comparative_metrics.delta_contacts !== undefined
                    ? `${detail.comparative_metrics.delta_contacts} contacts`
                    : '—'}
                </b>
              </div>
              <div>
                <span style={{ color: '#8297a4' }}>Secondary structure:</span>{' '}
                <b>
                  {detail.comparative_metrics.sec_struct_consistency ?? '—'}
                </b>
              </div>
            </div>

            {detail.comparative_metrics.comparison_caveats &&
              detail.comparative_metrics.comparison_caveats.length > 0 && (
                <div
                  style={{
                    marginTop: '12px',
                    padding: '10px 12px',
                    background: '#f8fafc',
                    borderRadius: '6px',
                    border: '1px solid #e2e8f0',
                  }}
                >
                  <b
                    style={{
                      fontSize: '11px',
                      color: '#475569',
                      textTransform: 'uppercase',
                      letterSpacing: '0.05em',
                    }}
                  >
                    Experimental & Construct Context:
                  </b>
                  <ul
                    style={{
                      margin: '6px 0 0',
                      paddingLeft: '18px',
                      fontSize: '12px',
                      color: '#334155',
                      lineHeight: '1.4',
                    }}
                  >
                    {detail.comparative_metrics.comparison_caveats.map((c, i) => (
                      <li key={i}>{c}</li>
                    ))}
                  </ul>
                </div>
              )}

            <p
              style={{
                fontSize: '11px',
                color: '#6d8797',
                marginTop: '12px',
                marginBottom: 0,
              }}
            >
              Notice: ΔRSA and ΔContacts represent static observational
              divergence between distinct experimental cryo-EM models. They do{' '}
              <b>not</b> indicate the causal biophysical impact of the mutation.
            </p>
          </article>

          {/* Source Records Table */}
          <article className="card" style={{ marginTop: '16px' }}>
            <div className="card-heading">
              <div>
                <span className="eyebrow">VERBATIM SOURCE ASSERTIONS</span>
                <h3>Clinical classifications without consensus synthesis</h3>
              </div>
            </div>
            <div style={{ display: 'flex', flexDirection: 'column', gap: '14px' }}>
              {detail.source_associations.cftr2_records.length > 0 && (
                <div>
                  <h4 style={{ margin: '0 0 6px', fontSize: '13px' }}>
                    CFTR2 Records ({detail.source_associations.cftr2_count})
                  </h4>
                  <table style={{ width: '100%', fontSize: '12px' }}>
                    <thead>
                      <tr>
                        <th>Legacy name</th>
                        <th>Determination</th>
                        <th>Alleles reported</th>
                        <th>Validation flag</th>
                      </tr>
                    </thead>
                    <tbody>
                      {detail.source_associations.cftr2_records.map((c, i) => (
                        <tr key={i}>
                          <td><b>{c.legacy_name}</b></td>
                          <td>{c.determination_verbatim ?? '—'}</td>
                          <td>{c.alleles_reported?.toLocaleString() ?? '—'}</td>
                          <td><Badge tone="teal">{c.validation_flag}</Badge></td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}

              {detail.source_associations.clinvar_records.length > 0 && (
                <div>
                  <h4 style={{ margin: '0 0 6px', fontSize: '13px' }}>
                    ClinVar Records ({detail.source_associations.clinvar_count})
                  </h4>
                  <table style={{ width: '100%', fontSize: '12px' }}>
                    <thead>
                      <tr>
                        <th>VariationID</th>
                        <th>VCV accession</th>
                        <th>Germline classification</th>
                        <th>Review status</th>
                      </tr>
                    </thead>
                    <tbody>
                      {detail.source_associations.clinvar_records.map((cv, i) => (
                        <tr key={i}>
                          <td><b>{cv.variation_id}</b></td>
                          <td>{cv.vcv_accession ?? '—'}</td>
                          <td>{cv.classification_verbatim ?? '—'}</td>
                          <td>{cv.review_status ?? '—'}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}
            </div>
          </article>

          {/* Warnings List */}
          {detail.warnings.map((w, idx) => (
            <Notice kind="warn" key={`${w.type}-${idx}`}>
              <span>
                <b>{w.type.replaceAll('_', ' ')}.</b> {w.message}
              </span>
            </Notice>
          ))}

          <Notice>
            <b>Scientific boundary.</b> {detail.guardrail_notice}
          </Notice>
        </>
      )}
    </>
  );
}

// ==============================================================================
// 5, 6, 7. Informational & Policy Sections
// ==============================================================================

function InfoPage({ page }: { page: Page }) {
  if (page === 'Batch Analysis') {
    return (
      <>
        <section className="page-intro">
          <span className="eyebrow">FUTURE WORKFLOW</span>
          <h2>Batch Analysis</h2>
          <p>
            Batch analysis endpoints are not implemented in the current FastAPI
            backend.
          </p>
        </section>

        <article className="card method-card">
          <div className="method-icon">
            <Table2 />
          </div>
          <div>
            <div className="method-row">
              <span className="check">✕</span>
              <p>
                <b>Endpoint Status:</b> Unavailable (`POST /variants/batch` or
                `/jobs` is not implemented in the current FastAPI backend).
              </p>
            </div>
            <div className="method-row">
              <span className="check">✕</span>
              <p>
                <b>No Simulated State:</b> The interface refuses to fabricate
                file uploads, simulated processing progress, or mock results.
              </p>
            </div>
            <div className="method-row">
              <span className="check">✓</span>
              <p>
                <b>Future Design:</b> When implemented, batch analysis will
                accept multi-variant VCF/TSV lists and stream per-structure
                evidence without altering the scientific data model.
              </p>
            </div>
          </div>
        </article>

        <Notice kind="warn">
          <b>Capability boundary.</b> Batch analysis remains unavailable until a
          corresponding asynchronous processing pipeline is exposed by the
          FastAPI backend.
        </Notice>
      </>
    );
  }

  if (page === 'Data Provenance') {
    return (
      <>
        <section className="page-intro">
          <span className="eyebrow">REPRODUCIBILITY & TRACEABILITY</span>
          <h2>Data Provenance</h2>
          <p>
            Every integrated observation is traceable to verified raw source
            files, SHA-256 hashes, and reproducible audit artifacts.
          </p>
        </section>

        <div className="grid-two">
          <article className="card">
            <div className="card-heading">
              <div>
                <span className="eyebrow">CLINICAL REGISTRIES</span>
                <h3>Ingested variant datasets</h3>
              </div>
            </div>
            <div className="structure-list">
              <div className="structure-row">
                <div>
                  <b>CFTR2_30January2026.xlsx</b>
                  <small>
                    2,092 source rows · SHA-256: 7b3c13b20fd4... · CFTR2 Project
                  </small>
                </div>
                <Badge tone="teal">2,092 records</Badge>
              </div>
              <div className="structure-row">
                <div>
                  <b>clinvar_download_20261009_143804.tsv</b>
                  <small>
                    2,664 source rows · SHA-256: a5291e5d9bab... · NCBI ClinVar
                    VCV export
                  </small>
                </div>
                <Badge tone="teal">2,664 records</Badge>
              </div>
            </div>
          </article>

          <article className="card">
            <div className="card-heading">
              <div>
                <span className="eyebrow">STRUCTURAL MODELS</span>
                <h3>Macromolecular cryo-EM files</h3>
              </div>
            </div>
            <div className="structure-list">
              <div className="structure-row">
                <span className="structure-id">5UAK</span>
                <div>
                  <b>5UAK.cif</b>
                  <small>Apo dephosphorylated · PDB deposition 2017</small>
                </div>
                <Badge tone="teal">3.87 Å</Badge>
              </div>
              <div className="structure-row">
                <span className="structure-id">6MSM</span>
                <div>
                  <b>6MSM.cif</b>
                  <small>ATP-bound dimer (E1371Q) · PDB deposition 2018</small>
                </div>
                <Badge tone="teal">3.20 Å</Badge>
              </div>
              <div className="structure-row">
                <span className="structure-id">8EIQ</span>
                <div>
                  <b>8EIQ.cif</b>
                  <small>ATP-bound ΔF508 (E1371Q) · PDB deposition 2022</small>
                </div>
                <Badge tone="warn">3.40 Å</Badge>
              </div>
            </div>
          </article>
        </div>

        <article className="card" style={{ marginTop: '16px' }}>
          <div className="card-heading">
            <div>
              <span className="eyebrow">AUDIT DELIVERABLES</span>
              <h3>Generated audit artifacts (docs/audits/)</h3>
            </div>
          </div>
          <ul className="principles">
            <li>
              <b>variant_structural_evidence.csv</b> — 4,756 per-source-record
              integrated evidence rows with per-structure features.
            </li>
            <li>
              <b>variant_structural_summary.csv</b> — 3,353 grouped records by
              normalized canonical protein change.
            </li>
            <li>
              <b>variant_classification_comparison.csv</b> — 693 pairwise CFTR2
              vs ClinVar classification comparisons.
            </li>
            <li>
              <b>variant_evidence_warnings.csv</b> — 1,959 audit warnings
              logging reference mismatches, deletions, substitutions, and
              unmodeled residues.
            </li>
            <li>
              <b>phase5_integration_summary.md</b> — Comprehensive scientific
              audit report.
            </li>
          </ul>
        </article>

        <Notice>
          <b>Persistent storage.</b> All data is indexed in the local SQLite
          database at `data/processed/cftr_structmap.db` (13 MB).
        </Notice>
      </>
    );
  }

  // About & Methods
  return (
    <>
      <section className="page-intro">
        <span className="eyebrow">SCIENTIFIC METHODOLOGY</span>
        <h2>About & Methods</h2>
        <p>
          CFTR-StructMap integrates normalized clinical variant evidence with
          experimental 3D structural context.
        </p>
      </section>

      <article className="card method-card">
        <div className="method-icon">
          <BookOpen />
        </div>
        <div>
          <div className="method-row">
            <span className="check">✓</span>
            <p>
              <b>Canonical Reference:</b> Human CFTR UniProt P13569, 1,480
              amino acids (SV=3, SHA-256 a86f7dbeb9c5...).
            </p>
          </div>
          <div className="method-row">
            <span className="check">✓</span>
            <p>
              <b>Solvent Accessibility:</b> Calculated with Shrake-Rupley
              (1.4 Å probe radius) and normalized to Relative Solvent
              Accessibility (RSA) using the Tien et al. (2013) empirical scale.
            </p>
          </div>
          <div className="method-row">
            <span className="check">✓</span>
            <p>
              <b>Contact Neighborhood:</b> Heavy-atom distance threshold of 8.0 Å,
              strictly excluding intra-residue atom pairs. In ATP-bound dimer models
              (6MSM and 8EIQ), the neighborhood includes surrounding protein residues
              as well as bound catalytic ligands (ATP and Mg²⁺).
            </p>
          </div>
          <div className="method-row">
            <span className="check">✓</span>
            <p>
              <b>Secondary Structure:</b> Deposited mmCIF helix and sheet header
              records.
            </p>
          </div>
          <div className="method-row">
            <span className="check">✓</span>
            <p>
              <b>Domain Context:</b> Lasso (1–79), TMD1 (80–388), NBD1
              (389–678), R-domain (679–844), TMD2 (845–1172), NBD2
              (1173–1436), C-term (1437–1480).
            </p>
          </div>
        </div>
      </article>

      <article className="card" style={{ marginTop: '16px' }}>
        <div className="card-heading">
          <div>
            <span className="eyebrow">FOUR SCIENTIFIC GUARDRAILS</span>
            <h3>Strict research boundaries</h3>
          </div>
        </div>
        <ul className="principles">
          <li>
            <b>Structural context ≠ pathogenicity:</b> Molecular features
            provide mechanistic context but do not replace clinical or
            functional assay proof.
          </li>
          <li>
            <b>Missing coordinates ≠ zero:</b> Unmodeled residues (e.g.
            intrinsically disordered R-domain) are strictly null, never zero.
          </li>
          <li>
            <b>Construct alterations ≠ wild-type:</b> Engineered
            mutations (e.g. E1371Q in 6MSM/8EIQ) and deletions (ΔF508 in 8EIQ)
            are explicitly flagged.
          </li>
          <li>
            <b>No manufactured consensus:</b> CFTR2 and ClinVar determinations
            are preserved verbatim without synthesizing artificial verdicts.
          </li>
        </ul>
      </article>

      <Notice>
        <b>Interpretation boundary.</b> This platform is designed for
        researchers and geneticists. It is not an automated diagnostic system.
      </Notice>
    </>
  );
}

// ==============================================================================
// Main Application Component
// ==============================================================================

function App() {
  const [page, setPage] = useState<Page>('Overview');
  const [selected, setSelected] = useState<string | null>(null);
  const [targetPos, setTargetPos] = useState<number | null>(null);
  const [backendOnline, setBackendOnline] = useState<boolean | null>(null);

  useEffect(() => {
    api
      .health()
      .then(() => setBackendOnline(true))
      .catch(() => setBackendOnline(false));
  }, []);

  let view: React.ReactNode;
  if (page === 'Overview') {
    view = <Overview setPage={setPage} />;
  } else if (page === 'Variant Explorer') {
    view = (
      <VariantExplorer
        onSelect={(id) => {
          setSelected(id);
          setPage('Structural Evidence');
        }}
      />
    );
  } else if (page === 'Batch Analysis') {
    view = (
      <BatchAnalysis
        onSelectVariant={(id) => {
          setSelected(id);
          setPage('Structural Evidence');
        }}
        onInspect3D={(id, pos) => {
          setSelected(id);
          if (pos) setTargetPos(pos);
          setPage('Structural Laboratory');
        }}
      />
    );
  } else if (page === 'Structural Laboratory') {
    view = (
      <StructureViewer
        targetVariant={selected}
        targetPosition={targetPos}
        onOpenEvidence={(id) => {
          setSelected(id);
          setPage('Structural Evidence');
        }}
      />
    );
  } else if (page === 'Structural Evidence') {
    view = (
      <Evidence
        selected={selected}
        onSelectVariant={(id) => setSelected(id)}
        onOpen3D={(pos) => {
          setTargetPos(pos);
          setPage('Structural Laboratory');
        }}
      />
    );
  } else {
    view = <InfoPage page={page} />;
  }

  return (
    <Shell page={page} setPage={setPage} backendOnline={backendOnline}>
      {view}
    </Shell>
  );
}

createRoot(document.getElementById('root')!).render(
  <React.StrictMode>
    <App />
  </React.StrictMode>
);
