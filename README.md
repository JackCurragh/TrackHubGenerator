# TrackHubGenerator

A Nextflow pipeline for converting BED12/Bedgraph files to  track hubs, originally designed for systematic processing of ribosome profiling data and translated region annotations.

## Overview

This pipeline automatically:
- Converts BED12 files to BigBed format (for translated regions)
- Converts Bedgraph files to BigWig format (for ribosome profiling signal)
- Creates organized trackhub directories for Ensembl FTP hosting
- Validates input data against formatting requirements
- Generates UCSC-compatible track hubs with proper metadata

## Directory Structure

```
trackhub-builder/
├── README.md                    # Main documentation
├── LICENSE                      # License file
├── main.nf                      # Main Nextflow workflow file
├── nextflow.config              # Nextflow configuration
├── modules/                     # Nextflow modules
│   ├── local/                   # Local modules
│   │   ├── bed_to_bigbed.nf     # BED to BigBed conversion 
│   │   ├── bedgraph_to_bigwig.nf # Bedgraph to BigWig conversion
│   │   └── create_trackhub.nf   # Track hub creation
│   └── nf-core/                 # Imported nf-core modules (if used)
├── conf/                        # Configuration directory
│   ├── base.config              # Base configuration
│   ├── test.config              # Test configuration
│   └── igenomes.config          # Reference genome configuration
├── bin/                         # Executable scripts directory
│   ├── validate_bed.py          # Validate BED12 format
│   ├── validate_bedgraph.py     # Validate Bedgraph format
│   ├── create_trackhub.py       # Python script to create trackhub
│   └── upload_trackhub.py       # Script to upload to Ensembl FTP
├── assets/                      # Static assets
│   ├── schema_input.json        # JSON schema for input validation
│   └── ensembl_colors.txt       # Standard Ensembl color schemes
├── docs/                        # Documentation
│   ├── usage.md                 # Usage documentation
│   ├── output.md                # Output documentation
│   └── examples.md              # Example commands
├── tests/                       # Tests directory for nf-test
│   ├── test_bed_conversion.nf.test # Test for BED conversion
│   ├── test_bedgraph_conversion.nf.test # Test for Bedgraph conversion
│   └── test_trackhub_creation.nf.test # Test for trackhub creation
├── .nf-test.yaml                # nf-test configuration
└── test_data/                   # Test data
    ├── bed/                     # Example BED12 files
    │   └── test_translons.bed   # Example translated regions
    └── bedgraph/                # Example Bedgraph files
        └── test_riboseq.bedgraph # Example ribosome profiling data
```

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/trackhub-builder.git
cd trackhub-builder

# Install Nextflow (if not already installed)
curl -fsSL get.nextflow.io | bash

# Move nextflow to a directory in your $PATH
mv nextflow ~/bin/
```

## Dependencies

- Nextflow (21.10.0 or higher)
- Python 3.8+
- Required Python packages:
  - `trackhub`
  - `pybedtools`
  - `pandas`
- UCSC tools:
  - `bedToBigBed`
  - `bedGraphToBigWig`
  - `fetchChromSizes`

You can install the Python dependencies with:

```bash
pip install trackhub pybedtools pandas
```

UCSC tools can be installed with:

```bash
# Download and install UCSC tools
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
chmod +x bedToBigBed bedGraphToBigWig fetchChromSizes
sudo mv bedToBigBed bedGraphToBigWig fetchChromSizes /usr/local/bin/
```

## Usage

### Basic Usage

```bash
nextflow run main.nf --input_dir /path/to/data --genome hg38 --track_type translons
```

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input_dir` | Directory containing input files | (required) |
| `--genome` | Reference genome (e.g., hg38, mm10) | (required) |
| `--track_type` | Type of tracks: 'translons' or 'riboseq' | (required) |
| `--phase` | Project phase (e.g., pilot, phase1) | 'pilot' |
| `--output_dir` | Output directory | './results' |
| `--hub_name` | Name for the trackhub | 'EnsemblTrackhub' |
| `--hub_short_label` | Short label for hub | 'Ensembl Trackhub' |
| `--hub_long_label` | Long description for hub | 'Ensembl Trackhub for Ribosome Profiling' |
| `--email` | Contact email | (required) |
| `--chrom_sizes` | Chromosome sizes file (if not using iGenomes) | null |
| `--upload` | Upload to Ensembl FTP when complete | false |
| `--ftp_user` | FTP username (if uploading) | null |
| `--ftp_pass` | FTP password (if uploading) | null |
| `--ftp_host` | FTP host (if uploading) | null |
| `--ftp_path` | FTP path (if uploading) | null |

## Input Files

### For `--track_type translons`:
- BED12 files containing translated regions (`.bed` or `.bed12` extension)

### For `--track_type riboseq`:
- Bedgraph files containing ribosome profiling signal (`.bedgraph` extension)

## Output Structure

The pipeline creates an organized directory structure ready for Ensembl FTP hosting:

```
output_dir/
├── translons/                # For translated regions
│   └── pilot/                # Or other phase
│       └── hg38/             # Reference genome
│           ├── hub.txt
│           ├── genomes.txt
│           ├── trackDb.txt
│           └── *.bb          # BigBed files
└── riboseq/                  # For ribosome profiling signal
    └── pilot/                # Or other phase
        └── hg38/             # Reference genome
            ├── hub.txt
            ├── genomes.txt
            ├── trackDb.txt
            └── *.bw          # BigWig files
```

## Examples

### Processing Translated Regions (BED12 files)

```bash
nextflow run main.nf \
  --input_dir /path/to/translon_data \
  --genome hg38 \
  --track_type translons \
  --phase pilot \
  --hub_name "TranslatedRegions" \
  --email "your.email@example.com"
```

### Processing Ribosome Profiling Signal (Bedgraph files)

```bash
nextflow run main.nf \
  --input_dir /path/to/riboseq_data \
  --genome hg38 \
  --track_type riboseq \
  --phase pilot \
  --hub_name "RibosomeProfileSignal" \
  --email "your.email@example.com"
```

### Uploading to Ensembl FTP Server

```bash
nextflow run main.nf \
  --input_dir /path/to/data \
  --genome hg38 \
  --track_type translons \
  --upload true \
  --ftp_user "username" \
  --ftp_pass "password" \
  --ftp_host "ftp.ensembl.org" \
  --ftp_path "/pub/trackhubs/translons/pilot"
```

## Testing

The pipeline includes tests using nf-test:

```bash
# Install nf-test if not already installed
wget -qO- https://code.askimed.com/install/nf-test | bash

# Run tests
nf-test test
```

## Schema Validation

The pipeline uses nf-schema for parameter validation:

```bash
# Install nf-schema if not already installed
pip install nf-schema

# Validate your parameters
nf-schema validate params.json
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the [MIT License](LICENSE).


### Destination FTP structure idea
```
ensembl_ftp/
├── translons/         # Contains BED12 files showing translated regions
│   ├── pilot/         # Initial data from consortium pilot phase
│   │   ├── hg38/      # Organized by reference genome
│   │   │   ├── hub.txt
│   │   │   ├── genomes.txt
│   │   │   └── trackDb.txt
│   │   └── mm10/
│   └── phase1/        # Future data releases
│
└── riboseq/           # Contains BigWig signal data from ribosome profiling
    ├── pilot/
    │   ├── hg38/
    │   │   ├── hub.txt
    │   │   ├── genomes.txt
    │   │   └── trackDb.txt
    │   └── mm10/
    └── phase1/
```