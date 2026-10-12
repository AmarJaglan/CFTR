import React, { useEffect, useRef, useState } from 'react';
import {
  AlertTriangle,
  Camera,
  CheckCircle2,
  Compass,
  ExternalLink,
  Layers,
  Maximize2,
  RotateCcw,
  Search,
  Sparkles,
} from 'lucide-react';
import { api, StructureMetadata, StructureResidueMapping } from '../lib/api';

declare global {
  interface Window {
    molstar?: {
      Viewer: {
        create: (
          element: HTMLElement | string,
          options?: Record<string, any>
        ) => Promise<any>;
      };
    };
  }
}

type StructureId = '5UAK' | '6MSM' | '8EIQ';

const STRUCTURE_INFO: Record<
  StructureId,
  {
    name: string;
    state: string;
    resolution: string;
    description: string;
    depositYear: string;
    ligands: string;
  }
> = {
  '5UAK': {
    name: '5UAK.cif',
    state: 'Apo, dephosphorylated monomer',
    resolution: '3.87 Å',
    description: 'Inward-facing, nucleotide-free reference state of human CFTR.',
    depositYear: '2017',
    ligands: 'No nucleotide bound',
  },
  '6MSM': {
    name: '6MSM.cif',
    state: 'ATP-bound dimer (E1371Q catalytic trap)',
    resolution: '3.20 Å',
    description: 'Active phosphorylated ATP-bound dimer trapped with catalytic E1371Q mutation.',
    depositYear: '2018',
    ligands: 'Bound ATP (ATP1504) and Magnesium (MG1502)',
  },
  '8EIQ': {
    name: '8EIQ.cif',
    state: 'ATP-bound dimer (ΔF508, E1371Q)',
    resolution: '3.40 Å',
    description: 'Active ATP-bound dimer carrying authentic patient mutation ΔF508.',
    depositYear: '2022',
    ligands: 'Bound ATP (ATP1504) and Magnesium (MG1502)',
  },
};

export function StructureViewer({
  initialStructure = '5UAK',
  targetVariant,
  targetPosition,
  onOpenEvidence,
}: {
  initialStructure?: StructureId;
  targetVariant?: string | null;
  targetPosition?: number | null;
  onOpenEvidence?: (variantId: string) => void;
}) {
  const [selectedStructure, setSelectedStructure] = useState<StructureId>(initialStructure);
  const [canonicalPos, setCanonicalPos] = useState<number>(targetPosition ?? 551);
  const [inputPosText, setInputPosText] = useState<string>(String(targetPosition ?? 551));

  // Viewer state
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerInstanceRef = useRef<any>(null);
  const [viewerLoading, setViewerLoading] = useState(true);
  const [viewerError, setViewerError] = useState<string | null>(null);

  // Mapping & Metadata
  const [mapping, setMapping] = useState<StructureResidueMapping | null>(null);
  const [metadata, setMetadata] = useState<StructureMetadata | null>(null);
  const [mappingLoading, setMappingLoading] = useState(false);

  // Sync target variant/position props if they change
  useEffect(() => {
    if (targetPosition) {
      setCanonicalPos(targetPosition);
      setInputPosText(String(targetPosition));
    }
  }, [targetPosition]);

  // Load structure metadata
  useEffect(() => {
    api
      .getStructureMetadata(selectedStructure)
      .then(setMetadata)
      .catch((err) => console.warn('Could not load structure metadata:', err));
  }, [selectedStructure]);

  // Load residue mapping whenever structure or position changes
  useEffect(() => {
    if (!canonicalPos || canonicalPos < 1 || canonicalPos > 1480) {
      setMapping(null);
      return;
    }

    setMappingLoading(true);
    api
      .getStructureResidue(selectedStructure, canonicalPos)
      .then((m) => setMapping(m))
      .catch(() => setMapping(null))
      .finally(() => setMappingLoading(false));
  }, [selectedStructure, canonicalPos]);

  // Initialize Mol* viewer instance on selected structure
  useEffect(() => {
    let isCancelled = false;

    // Clean up previous viewer if exists
    if (viewerInstanceRef.current) {
      try {
        viewerInstanceRef.current.dispose();
      } catch (err) {
        console.warn('Error disposing Mol* viewer:', err);
      }
      viewerInstanceRef.current = null;
    }

    if (!containerRef.current) return;
    containerRef.current.innerHTML = '';

    setViewerLoading(true);
    setViewerError(null);

    // Wait for window.molstar to be ready
    const initViewer = async () => {
      let attempts = 0;
      while (!window.molstar && attempts < 20) {
        await new Promise((r) => setTimeout(r, 100));
        attempts++;
      }

      if (!window.molstar) {
        if (!isCancelled) {
          setViewerError('Mol* library not loaded. Check network or bundle configuration.');
          setViewerLoading(false);
        }
        return;
      }

      try {
        const viewer = await window.molstar.Viewer.create(containerRef.current!, {
          layoutIsExpanded: false,
          layoutShowControls: false,
          layoutShowRemoteState: false,
          layoutShowSequence: true,
          layoutShowLog: false,
          layoutShowLeftPanel: false,
          viewportShowExpand: false,
          viewportShowSelectionMode: false,
          viewportShowAnimation: false,
          collapseLeftPanel: true,
          pdbProvider: 'pdbe',
          emdbProvider: 'pdbe',
        });

        if (isCancelled) {
          viewer.dispose();
          return;
        }

        viewerInstanceRef.current = viewer;

        // Load mmCIF file from secure API endpoint
        const coordUrl = api.getStructureCoordinatesUrl(selectedStructure);
        await viewer.loadStructureFromUrl(coordUrl, 'mmcif');

        if (!isCancelled) {
          setViewerLoading(false);
        }
      } catch (err: unknown) {
        if (!isCancelled) {
          const msg = err instanceof Error ? err.message : 'Error rendering 3D coordinates';
          setViewerError(msg);
          setViewerLoading(false);
        }
      }
    };

    initViewer();

    return () => {
      isCancelled = true;
      if (viewerInstanceRef.current) {
        try {
          viewerInstanceRef.current.dispose();
        } catch {
          // ignore cleanup errors
        }
        viewerInstanceRef.current = null;
      }
    };
  }, [selectedStructure]);

  // Camera reset
  const handleResetCamera = () => {
    if (!viewerInstanceRef.current) return;
    try {
      const plugin = viewerInstanceRef.current.plugin;
      if (plugin?.canvas3d) {
        plugin.canvas3d.requestCameraReset();
      }
    } catch (e) {
      console.warn('Camera reset error:', e);
    }
  };

  const handlePositionSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    const p = parseInt(inputPosText, 10);
    if (!isNaN(p) && p >= 1 && p <= 1480) {
      setCanonicalPos(p);
    }
  };

  const currentInfo = STRUCTURE_INFO[selectedStructure];

  return (
    <>
      <section className="page-intro">
        <span className="eyebrow">EXPERIMENTAL MACROMOLECULAR MODELS</span>
        <h2>Structural Laboratory & 3D Viewer</h2>
        <p>
          Inspect real experimental Cryo-EM coordinates for human CFTR (5UAK, 6MSM, 8EIQ),
          explore conformational changes, and map canonical variants onto 3D residue positions.
        </p>
      </section>

      {/* Structure Selector Bar */}
      <article className="card" style={{ marginBottom: '16px', padding: '12px 16px' }}>
        <div
          style={{
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'center',
            flexWrap: 'wrap',
            gap: '12px',
          }}
        >
          <div style={{ display: 'flex', gap: '8px', alignItems: 'center' }}>
            <span style={{ fontSize: '13px', fontWeight: 600, color: '#334155' }}>
              Select Model:
            </span>
            {(['5UAK', '6MSM', '8EIQ'] as StructureId[]).map((id) => (
              <button
                key={id}
                className={selectedStructure === id ? 'primary' : 'secondary'}
                style={{
                  padding: '6px 14px',
                  fontSize: '13px',
                  fontWeight: selectedStructure === id ? 700 : 500,
                }}
                onClick={() => setSelectedStructure(id)}
              >
                <b>{id}</b> ({STRUCTURE_INFO[id].resolution})
              </button>
            ))}
          </div>

          <div style={{ display: 'flex', gap: '8px', alignItems: 'center' }}>
            <button
              className="secondary"
              style={{ padding: '6px 12px', fontSize: '12px' }}
              onClick={handleResetCamera}
              title="Reset 3D camera to default orientation"
            >
              <RotateCcw style={{ width: '13px', marginRight: '4px', display: 'inline-block' }} />
              Reset Camera
            </button>
          </div>
        </div>
      </article>

      {/* Main 2-Column Workspace */}
      <div
        style={{
          display: 'grid',
          gridTemplateColumns: 'minmax(0, 1.8fr) minmax(320px, 1fr)',
          gap: '16px',
          alignItems: 'start',
        }}
      >
        {/* 3D Canvas Card */}
        <article
          className="card"
          style={{
            padding: 0,
            overflow: 'hidden',
            position: 'relative',
            borderRadius: '8px',
            background: '#0f172a',
          }}
        >
          {/* Status Overlay */}
          {viewerLoading && (
            <div
              style={{
                position: 'absolute',
                top: 0,
                left: 0,
                right: 0,
                bottom: 0,
                background: 'rgba(15, 23, 42, 0.85)',
                color: '#f8fafc',
                display: 'flex',
                flexDirection: 'column',
                justifyContent: 'center',
                alignItems: 'center',
                zIndex: 10,
              }}
            >
              <Compass
                style={{
                  width: '36px',
                  height: '36px',
                  color: '#14b8a6',
                  animation: 'spin 2s linear infinite',
                  marginBottom: '12px',
                }}
              />
              <p style={{ margin: 0, fontWeight: 600 }}>Streaming {selectedStructure} Coordinates...</p>
              <small style={{ color: '#94a3b8', marginTop: '4px' }}>
                Parsing mmCIF atoms and initializing WebGL shaders
              </small>
            </div>
          )}

          {viewerError && (
            <div
              style={{
                position: 'absolute',
                top: 0,
                left: 0,
                right: 0,
                bottom: 0,
                background: '#fef2f2',
                color: '#991b1b',
                display: 'flex',
                flexDirection: 'column',
                justifyContent: 'center',
                alignItems: 'center',
                padding: '24px',
                zIndex: 10,
              }}
            >
              <AlertTriangle style={{ width: '32px', height: '32px', marginBottom: '8px' }} />
              <p style={{ fontWeight: 600, margin: '0 0 4px' }}>Viewer Initialization Notice</p>
              <p style={{ fontSize: '13px', textAlign: 'center', margin: 0 }}>{viewerError}</p>
            </div>
          )}

          {/* Interactive Mol* Canvas Container */}
          <div
            ref={containerRef}
            style={{
              width: '100%',
              height: '540px',
              position: 'relative',
            }}
          />

          {/* Canvas Footer Bar */}
          <div
            style={{
              padding: '8px 12px',
              background: '#1e293b',
              color: '#94a3b8',
              fontSize: '11px',
              display: 'flex',
              justifyContent: 'space-between',
              alignItems: 'center',
            }}
          >
            <span>
              Controls: Left-click rotate · Right-click pan · Scroll zoom · Sequence bar top
            </span>
            <span style={{ color: '#38bdf8' }}>
              {selectedStructure} ({currentInfo.resolution})
            </span>
          </div>
        </article>

        {/* Right Inspection & Residue Mapping Panel */}
        <div>
          {/* Target Residue Selector */}
          <article className="card" style={{ marginBottom: '16px' }}>
            <h3 style={{ margin: '0 0 10px', fontSize: '15px' }}>Residue Inspection</h3>

            {/* Quick Benchmark Buttons */}
            <div style={{ display: 'flex', gap: '6px', flexWrap: 'wrap', marginBottom: '12px' }}>
              {[
                { label: 'G551 (NBD1)', pos: 551 },
                { label: 'F508 (NBD1)', pos: 508 },
                { label: 'E1371 (NBD2)', pos: 1371 },
                { label: 'M470 (NBD1)', pos: 470 },
                { label: 'R117 (TMD1)', pos: 117 },
              ].map((b) => (
                <button
                  key={b.pos}
                  className="secondary"
                  style={{
                    padding: '4px 8px',
                    fontSize: '11px',
                    background: canonicalPos === b.pos ? '#e2e8f0' : undefined,
                    fontWeight: canonicalPos === b.pos ? 700 : 400,
                  }}
                  onClick={() => {
                    setCanonicalPos(b.pos);
                    setInputPosText(String(b.pos));
                  }}
                >
                  {b.label}
                </button>
              ))}
            </div>

            {/* Position Search Form */}
            <form onSubmit={handlePositionSubmit} style={{ display: 'flex', gap: '8px' }}>
              <input
                type="number"
                min="1"
                max="1480"
                value={inputPosText}
                onChange={(e) => setInputPosText(e.target.value)}
                placeholder="Canonical residue 1–1480"
                style={{
                  flex: 1,
                  padding: '6px 10px',
                  fontSize: '13px',
                  borderRadius: '6px',
                  border: '1px solid #cbd5e1',
                }}
              />
              <button className="primary" type="submit" style={{ padding: '6px 14px', fontSize: '13px' }}>
                Map Position
              </button>
            </form>
          </article>

          {/* Residue Mapping & Guardrails Card */}
          <article className="card" style={{ marginBottom: '16px' }}>
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '10px' }}>
              <span className="eyebrow">CANONICAL CORRESPONDENCE</span>
              <span
                className={`badge ${
                  mapping?.mapping_status === 'MAPPED_OBSERVED'
                    ? 'good'
                    : mapping?.mapping_status === 'NOT_COVERED'
                      ? 'warn'
                      : 'neutral'
                }`}
              >
                {mapping?.mapping_status ?? 'NOT_MAPPED'}
              </span>
            </div>

            <h4 style={{ margin: '0 0 12px', fontSize: '16px' }}>
              Position {canonicalPos} ({mapping?.canonical_amino_acid ?? '—'})
            </h4>

            {mapping ? (
              <dl style={{ margin: 0 }}>
                <dt>Structure chain</dt>
                <dd>Chain {mapping.auth_asym_id} ({mapping.label_asym_id})</dd>

                <dt>Author residue number</dt>
                <dd>{mapping.auth_seq_id ?? 'None'}</dd>

                <dt>Observed amino acid</dt>
                <dd>
                  {mapping.structure_amino_acid ? (
                    <b>{mapping.structure_amino_acid}</b>
                  ) : (
                    <span style={{ color: '#94a3b8' }}>— (Unmodeled)</span>
                  )}
                </dd>

                <dt>Cα atom present</dt>
                <dd>{mapping.has_ca_atom ? '✓ Yes' : '✕ No'}</dd>
              </dl>
            ) : (
              <p style={{ color: '#64748b', fontSize: '13px' }}>
                {mappingLoading ? 'Loading residue mapping...' : 'No mapping recorded for this position.'}
              </p>
            )}

            {/* Scientific Guardrails & Specific Contexts */}
            {canonicalPos === 551 && (
              <div
                style={{
                  marginTop: '12px',
                  padding: '8px 10px',
                  background: '#f0fdfa',
                  borderRadius: '6px',
                  border: '1px solid #ccfbf1',
                  fontSize: '11px',
                  color: '#0f766e',
                }}
              >
                <b>Wild-Type Reference State:</b> The Cryo-EM models contain reference Glycine 551 (LSGGQ motif).
                When evaluating G551D, the viewer displays the wild-type Glycine context, not an in silico mutant sidechain.
              </div>
            )}

            {canonicalPos === 508 && selectedStructure === '8EIQ' && (
              <div
                style={{
                  marginTop: '12px',
                  padding: '8px 10px',
                  background: '#fffbeb',
                  borderRadius: '6px',
                  border: '1px solid #fef3c7',
                  fontSize: '11px',
                  color: '#92400e',
                }}
              >
                <b>Construct Deletion:</b> Residue F508 is deleted from the 8EIQ construct (authentic ΔF508 patient mutation).
                Atomic coordinates are absent. Switch to <b>5UAK</b> or <b>6MSM</b> to observe wild-type Phe 508.
              </div>
            )}

            {canonicalPos === 1371 && ['6MSM', '8EIQ'].includes(selectedStructure) && (
              <div
                style={{
                  marginTop: '12px',
                  padding: '8px 10px',
                  background: '#eff6ff',
                  borderRadius: '6px',
                  border: '1px solid #dbeafe',
                  fontSize: '11px',
                  color: '#1e40af',
                }}
              >
                <b>Engineered Mutation:</b> Models 6MSM and 8EIQ carry the catalytic substitution E1371Q introduced
                to trap the ATP-bound dimer by inhibiting ATP hydrolysis.
              </div>
            )}

            {onOpenEvidence && (
              <button
                className="secondary"
                style={{ marginTop: '14px', width: '100%', fontSize: '12px', padding: '6px' }}
                onClick={() => onOpenEvidence(String(canonicalPos))}
              >
                <ExternalLink style={{ width: '12px', marginRight: '4px', display: 'inline-block' }} />
                Open in Structural Evidence Table
              </button>
            )}
          </article>

          {/* Experimental Model Metadata Card */}
          <article className="card">
            <span className="eyebrow">CRYO-EM DEPOSITION CONTEXT</span>
            <h4 style={{ margin: '4px 0 8px', fontSize: '14px' }}>{currentInfo.name}</h4>
            <p style={{ margin: '0 0 10px', fontSize: '12px', color: '#475569' }}>
              {currentInfo.description}
            </p>
            <dl style={{ margin: 0, fontSize: '12px' }}>
              <dt>Conformational state</dt>
              <dd>{currentInfo.state}</dd>
              <dt>Nominal resolution</dt>
              <dd><b>{currentInfo.resolution}</b></dd>
              <dt>Deposition year</dt>
              <dd>{currentInfo.depositYear}</dd>
              <dt>Ligand context</dt>
              <dd>{currentInfo.ligands}</dd>
            </dl>
          </article>
        </div>
      </div>
    </>
  );
}
