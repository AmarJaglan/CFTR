import React from 'react';
import { describe, expect, it, vi } from 'vitest';
import { render, screen, fireEvent, waitFor } from '@testing-library/react';
import { BatchAnalysis } from './BatchAnalysis';
import { api, BatchVariantResponse } from '../lib/api';

describe('BatchAnalysis component', () => {
  it('renders initial input interface and examples button', () => {
    const onSelect = vi.fn();
    const onInspect = vi.fn();
    render(<BatchAnalysis onSelectVariant={onSelect} onInspect3D={onInspect} />);

    expect(screen.getByText('Batch Variant Analysis')).toBeDefined();
    expect(screen.getByText('Paste Identifiers')).toBeDefined();
    expect(screen.getByText('Upload CSV File')).toBeDefined();
    expect(screen.getByText('Load Test Examples')).toBeDefined();
  });

  it('prefills examples when Load Test Examples is clicked', () => {
    render(<BatchAnalysis onSelectVariant={vi.fn()} onInspect3D={vi.fn()} />);
    const exampleBtn = screen.getByText('Load Test Examples');
    fireEvent.click(exampleBtn);

    const textarea = screen.getByPlaceholderText(/Enter one variant per line/i) as HTMLTextAreaElement;
    expect(textarea.value).toContain('G551D');
    expect(textarea.value).toContain('F508del');
    expect(textarea.value).toContain('p.Gly551Asp');
  });

  it('submits batch variants and displays results table and summary metrics', async () => {
    const mockResponse: BatchVariantResponse = {
      batch_id: 'batch-test-001',
      timestamp: '2026-10-12T12:00:00Z',
      total_input_rows: 3,
      unique_identifiers: 3,
      resolved_count: 2,
      unresolved_count: 1,
      evidence_available_count: 2,
      results: [
        {
          row_number: 1,
          original_input: 'G551D',
          normalized_identifier: 'G551D',
          status: 'RESOLVED',
          canonical_protein: 'G551D',
          canonical_position: 551,
          ref_aa: 'G',
          alt_aa: 'D',
          domain_name: 'NBD1',
          source_record_count: 2,
          source_determinations: ['CF-causing'],
          has_structural_evidence: true,
          structures: {
            '5UAK': {
              state: 'apo',
              mapping_status: 'MAPPED_OBSERVED',
              structure_aa: 'G',
              feature_status: 'ELIGIBLE',
              failure_reason: null,
              sasa_angstrom2: 8.0,
              rsa_ratio: 0.0769,
              exposure_category: 'BURIED',
              total_contacts: 19,
              intra_contacts: 19,
              inter_contacts: 0,
              secondary_structure: 'HELIX',
              construct_note: null,
            },
            '6MSM': {
              state: 'dimer',
              mapping_status: 'MAPPED_OBSERVED',
              structure_aa: 'G',
              feature_status: 'ELIGIBLE',
              failure_reason: null,
              sasa_angstrom2: 0.0,
              rsa_ratio: 0.0,
              exposure_category: 'BURIED',
              total_contacts: 30,
              intra_contacts: 30,
              inter_contacts: 0,
              secondary_structure: 'HELIX',
              construct_note: null,
            },
            '8EIQ': {
              state: 'dimer',
              mapping_status: 'MAPPED_OBSERVED',
              structure_aa: 'G',
              feature_status: 'ELIGIBLE',
              failure_reason: null,
              sasa_angstrom2: 0.0,
              rsa_ratio: 0.0,
              exposure_category: 'BURIED',
              total_contacts: 30,
              intra_contacts: 30,
              inter_contacts: 0,
              secondary_structure: 'HELIX',
              construct_note: null,
            },
          },
          comparative_metrics: {
            observation_pattern: 'ALL_THREE_OBSERVED',
            construct_variation: 'CANONICAL_WT',
            is_fully_comparable: true,
            delta_rsa: 0.0769,
            delta_contacts: 11,
            sec_struct_consistency: 'IDENTICAL',
            comparison_notes: null,
          },
          warnings: [],
          failure_reason: null,
        },
        {
          row_number: 2,
          original_input: 'F508del',
          normalized_identifier: 'F508del',
          status: 'RESOLVED',
          canonical_protein: 'F508del',
          canonical_position: 508,
          ref_aa: 'F',
          alt_aa: 'del',
          domain_name: 'NBD1',
          source_record_count: 1,
          source_determinations: ['CF-causing'],
          has_structural_evidence: true,
          structures: {
            '5UAK': {
              state: 'apo',
              mapping_status: 'MAPPED_OBSERVED',
              structure_aa: 'F',
              feature_status: 'ELIGIBLE',
              failure_reason: null,
              sasa_angstrom2: 35.0,
              rsa_ratio: 0.16,
              exposure_category: 'BURIED',
              total_contacts: 24,
              intra_contacts: 24,
              inter_contacts: 0,
              secondary_structure: 'HELIX',
              construct_note: null,
            },
            '6MSM': {
              state: 'dimer',
              mapping_status: 'MAPPED_OBSERVED',
              structure_aa: 'F',
              feature_status: 'ELIGIBLE',
              failure_reason: null,
              sasa_angstrom2: 28.0,
              rsa_ratio: 0.13,
              exposure_category: 'BURIED',
              total_contacts: 26,
              intra_contacts: 26,
              inter_contacts: 0,
              secondary_structure: 'HELIX',
              construct_note: null,
            },
            '8EIQ': {
              state: 'dimer',
              mapping_status: 'NOT_COVERED',
              structure_aa: null,
              feature_status: 'INELIGIBLE_NOT_COVERED',
              failure_reason: 'CONSTRUCT_DELETION',
              sasa_angstrom2: null,
              rsa_ratio: null,
              exposure_category: null,
              total_contacts: null,
              intra_contacts: null,
              inter_contacts: null,
              secondary_structure: null,
              construct_note: 'Construct deletion: Delta-F508',
            },
          },
          comparative_metrics: null,
          warnings: [],
          failure_reason: null,
        },
        {
          row_number: 3,
          original_input: 'UnknownXYZ999',
          normalized_identifier: null,
          status: 'UNSUPPORTED_FORMAT',
          canonical_protein: null,
          canonical_position: null,
          ref_aa: null,
          alt_aa: null,
          domain_name: null,
          source_record_count: 0,
          source_determinations: [],
          has_structural_evidence: false,
          structures: null,
          comparative_metrics: null,
          warnings: [],
          failure_reason: "Unable to resolve 'UnknownXYZ999' as a supported CFTR variant identifier.",
        },
      ],
    };

    const batchSpy = vi.spyOn(api, 'batchVariants').mockResolvedValue(mockResponse);
    const onSelect = vi.fn();
    const onInspect = vi.fn();

    render(<BatchAnalysis onSelectVariant={onSelect} onInspect3D={onInspect} />);

    // Load examples and submit
    fireEvent.click(screen.getByText('Load Test Examples'));
    fireEvent.click(screen.getByRole('button', { name: /Run Batch Analysis/i }));

    await waitFor(() => {
      expect(screen.getByText('Batch Results')).toBeDefined();
    });

    // Check summary metrics rendered
    expect(screen.getByText('Input Rows')).toBeDefined();
    expect(screen.getByText('Found in database')).toBeDefined();

    // Check results table rows
    expect(screen.getByText('UnknownXYZ999')).toBeDefined();
    expect(screen.getByText('UNSUPPORTED_FORMAT')).toBeDefined();
    expect(screen.getByText('Export CSV')).toBeDefined();

    batchSpy.mockRestore();
  });
});
