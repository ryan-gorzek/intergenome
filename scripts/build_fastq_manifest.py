#!/usr/bin/env python3
"""
build_fastq_manifest.py

Convert NeMO archive manifest TSVs into pipeline-ready format for DOWNLOAD_FASTQ process.

Usage:
    python build_fastq_manifest.py [OPTIONS] MANIFEST_FILES...

Arguments:
    MANIFEST_FILES  - One or more NeMO manifest TSV files (from fetch_nemo_manifests.sh)

Options:
    -o, --output FILE    - Output manifest file (default: manifests/fastq_manifest.tsv)
    -n, --max-samples N  - Limit to N samples per species (default: all)
    -s, --species LIST   - Filter to specific species (comma-separated)
    --validate           - Check that URLs are accessible (slow)
    -h, --help           - Show this help message

Output format (TSV):
    URL<TAB>folder<TAB>sample

Example:
    python build_fastq_manifest.py manifests/nemo/*.tsv -o manifests/fastq_manifest.tsv
    python build_fastq_manifest.py manifests/nemo/human_M1_manifest.tsv -n 5
"""

import argparse
import sys
from pathlib import Path
from urllib.request import urlopen
from urllib.error import URLError


def parse_nemo_manifest(filepath: Path) -> list[dict]:
    """Parse a NeMO manifest TSV file and extract FASTQ entries."""
    entries = []

    with open(filepath, 'r') as f:
        header = f.readline().strip().split('\t')

        # Find relevant columns
        col_indices = {}
        for i, col in enumerate(header):
            col_lower = col.lower().replace(' ', '_')
            if 'file_name' in col_lower or col_lower == 'filename':
                col_indices['filename'] = i
            elif 'released_url' in col_lower or col_lower == 'url':
                col_indices['url'] = i
            elif 'sample_id' in col_lower or col_lower == 'sample':
                col_indices['sample'] = i
            elif col_lower == 'species':
                col_indices['species'] = i

        # Also try common column names directly
        for i, col in enumerate(header):
            if col == 'File_name':
                col_indices['filename'] = i
            elif col == 'Released_URL':
                col_indices['url'] = i
            elif col == 'Sample_ID':
                col_indices['sample'] = i
            elif col == 'Species':
                col_indices['species'] = i

        if 'url' not in col_indices:
            print(f"[WARN] No URL column found in {filepath}", file=sys.stderr)
            print(f"       Columns: {header}", file=sys.stderr)
            return entries

        for line in f:
            fields = line.strip().split('\t')
            if len(fields) <= col_indices['url']:
                continue

            url = fields[col_indices['url']]
            if not url or not url.startswith('http'):
                continue

            # Only include FASTQ files
            if '.fastq' not in url.lower():
                continue

            entry = {'url': url}

            if 'filename' in col_indices and len(fields) > col_indices['filename']:
                entry['filename'] = fields[col_indices['filename']]
            else:
                entry['filename'] = url.split('/')[-1]

            if 'sample' in col_indices and len(fields) > col_indices['sample']:
                entry['sample'] = fields[col_indices['sample']]
            else:
                # Extract sample from filename
                entry['sample'] = entry['filename'].replace('.fastq.tar', '').replace('.fastq.gz', '')

            if 'species' in col_indices and len(fields) > col_indices['species']:
                entry['species'] = fields[col_indices['species']]
            else:
                # Try to infer from filepath
                entry['species'] = filepath.stem.replace('_M1_manifest', '')

            entries.append(entry)

    return entries


def validate_url(url: str) -> bool:
    """Check if URL is accessible (HEAD request)."""
    try:
        req = urlopen(url, timeout=10)
        req.close()
        return True
    except URLError:
        return False


def main():
    parser = argparse.ArgumentParser(
        description='Convert NeMO manifests to pipeline-ready format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('manifests', nargs='+', type=Path,
                        help='NeMO manifest TSV files')
    parser.add_argument('-o', '--output', type=Path,
                        default=Path('manifests/fastq_manifest.tsv'),
                        help='Output manifest file')
    parser.add_argument('-n', '--max-samples', type=int, default=None,
                        help='Limit samples per species')
    parser.add_argument('-s', '--species', type=str, default=None,
                        help='Filter to specific species (comma-separated)')
    parser.add_argument('--validate', action='store_true',
                        help='Validate URLs are accessible')

    args = parser.parse_args()

    # Parse species filter
    species_filter = None
    if args.species:
        species_filter = set(s.strip().lower() for s in args.species.split(','))

    # Collect all entries
    all_entries = []
    for manifest_path in args.manifests:
        if not manifest_path.exists():
            print(f"[WARN] File not found: {manifest_path}", file=sys.stderr)
            continue

        print(f"[INFO] Parsing: {manifest_path}")
        entries = parse_nemo_manifest(manifest_path)
        print(f"       Found {len(entries)} FASTQ entries")
        all_entries.extend(entries)

    if not all_entries:
        print("[ERROR] No entries found in manifests", file=sys.stderr)
        sys.exit(1)

    # Filter by species
    if species_filter:
        all_entries = [e for e in all_entries
                       if e.get('species', '').lower() in species_filter]
        print(f"[INFO] After species filter: {len(all_entries)} entries")

    # Limit samples per species
    if args.max_samples:
        by_species = {}
        for e in all_entries:
            sp = e.get('species', 'unknown')
            if sp not in by_species:
                by_species[sp] = []
            by_species[sp].append(e)

        all_entries = []
        for sp, entries in by_species.items():
            all_entries.extend(entries[:args.max_samples])
            if len(entries) > args.max_samples:
                print(f"[INFO] Limited {sp} to {args.max_samples} samples (had {len(entries)})")

    # Validate URLs if requested
    if args.validate:
        print("[INFO] Validating URLs (this may take a while)...")
        valid_entries = []
        for e in all_entries:
            if validate_url(e['url']):
                valid_entries.append(e)
            else:
                print(f"[WARN] Invalid URL: {e['url']}", file=sys.stderr)
        all_entries = valid_entries
        print(f"[INFO] {len(all_entries)} valid URLs")

    # Write output manifest
    args.output.parent.mkdir(parents=True, exist_ok=True)

    with open(args.output, 'w') as f:
        for e in all_entries:
            # Format: URL<TAB>folder<TAB>sample
            # folder is species-based for organization
            folder = f"U01_Lein/{e.get('species', 'unknown')}"
            sample = e['sample']
            url = e['url']
            f.write(f"{url}\t{folder}\n")

    print(f"[INFO] Wrote {len(all_entries)} entries to: {args.output}")
    print(f"")
    print(f"Next step: Run nextflow with --fastq_mode download --fastq_manifest {args.output}")


if __name__ == '__main__':
    main()
