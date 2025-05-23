#!/usr/bin/env python3
"""
Create UCSC and Ensembl compatible track hubs from globbed files with sample metadata.

This script generates track hubs for genomic data files (BigWig or BigBed) by matching
glob patterns against a sample sheet that defines metadata and display options.
"""

import os
import json
import csv
import glob
import re
import argparse
import logging
from pathlib import Path
from typing import Dict, List, Optional, Pattern, Any
import trackhub

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


# Dictionary of track parameters that work in both UCSC and Ensembl
COMPATIBLE_TRACK_PARAMS = {
    'bigWig': {
        'common': ['visibility', 'color', 'priority', 'autoScale', 'alwaysZero', 'maxHeightPixels'],
        'ucsc_only': ['graphType', 'viewLimits', 'smoothingWindow'],
        'ensembl_only': ['type', 'format']
    },
    'bigBed': {
        'common': ['visibility', 'color', 'priority'],
        'ucsc_only': ['scoreFilter', 'itemRgb', 'spectrum'],
        'ensembl_only': ['type', 'format']
    }
}


class SampleMetadata:
    """Represents sample metadata from the sample sheet."""
    
    def __init__(self, row: Dict[str, str]):
        """
        Initialize sample metadata from a sample sheet row.
        
        Args:
            row: Dictionary representing a row from the sample sheet
        """
        self.sample_id: str = row['sample_id']
        
        # Optional fields with defaults
        self.short_label: str = row.get('short_label', self.sample_id)
        self.long_label: str = row.get('long_label', self.sample_id)
        self.color: str = row.get('color', '')
        self.visibility: str = row.get('visibility', '')
        self.priority: str = row.get('priority', '1')
        
        # Additional parameters as JSON
        self.additional_params: Dict[str, Any] = {}
        if 'additional_params' in row and row['additional_params']:
            try:
                self.additional_params = json.loads(row['additional_params'])
            except json.JSONDecodeError:
                logger.warning(f"Could not parse additional_params for {self.sample_id}: {row['additional_params']}")


class TrackFile:
    """Represents a discovered track file with extracted metadata."""
    
    def __init__(
        self, 
        file_path: str, 
        track_type: str,
        sample_pattern: Optional[Pattern] = None,
        annotation_pattern: Optional[Pattern] = None
    ):
        """
        Initialize a track file with extracted metadata.
        
        Args:
            file_path: Path to the track file
            track_type: Type of track (bigwig or bigbed)
            sample_pattern: Regex pattern with named group 'sample_id' to extract sample ID
            annotation_pattern: Regex pattern with named group 'annotation_type' to extract annotation type
        """
        self.file_path: str = file_path
        self.track_type: str = track_type.lower()
        self.filename: str = os.path.basename(file_path)
        self.basename: str = os.path.splitext(self.filename)[0]
        
        # Extract sample ID from filename
        self.sample_id: str = self._extract_sample_id(sample_pattern)
        
        # Extract annotation type for bigbed files
        self.annotation_type: str = ''
        if self.track_type == 'bigbed':
            self.annotation_type = self._extract_annotation_type(annotation_pattern)
        
        # Track display properties will be set from sample metadata
        self.short_label: str = self.basename
        self.long_label: str = self.basename
        self.color: str = ''
        self.visibility: str = ''
        self.priority: str = '1'
        self.additional_params: Dict[str, Any] = {}
    
    def _extract_sample_id(self, pattern: Optional[Pattern] = None) -> str:
        """Extract sample ID from filename using pattern or default method."""
        if pattern:
            match = pattern.search(self.filename)
            if match and 'sample_id' in match.groupdict():
                return match.group('sample_id')
        
        # Default: use the part before the first dot or underscore
        return self.basename.split('.')[0].split('_')[0]
    
    def _extract_annotation_type(self, pattern: Optional[Pattern] = None) -> str:
        """Extract annotation type from filename using pattern or default method."""
        if pattern:
            match = pattern.search(self.filename)
            if match and 'annotation_type' in match.groupdict():
                return match.group('annotation_type')
        
        # Default: try to extract a tool name after the sample ID
        parts = self.basename.split('_')
        if len(parts) > 1:
            return parts[1]
        
        return 'default'
    
    def apply_metadata(self, metadata: SampleMetadata) -> None:
        """
        Apply sample metadata to this track file.
        
        Args:
            metadata: SampleMetadata object to apply
        """
        self.short_label = metadata.short_label
        self.long_label = metadata.long_label
        
        # Only apply color if not already set
        if not self.color and metadata.color:
            self.color = metadata.color
        
        # Only apply visibility if not already set
        if not self.visibility and metadata.visibility:
            self.visibility = metadata.visibility
        
        # Only apply priority if not already set
        if not self.priority and metadata.priority:
            self.priority = metadata.priority
        
        # Merge additional parameters
        self.additional_params.update(metadata.additional_params)
    
    def get_track_params(self, parent=None, ensembl_compatible: bool = False) -> Dict[str, Any]:
        """
        Get track parameters for creating a trackhub Track object.
        
        Args:
            parent: Optional parent track for composite/supertrack organization
            ensembl_compatible: If True, include parameters for Ensembl compatibility
            
        Returns:
            Dictionary of track parameters
        """
        # Generate unique track name
        track_name = f"{self.track_type}_{self.basename}"
        
        # Set default color if not specified
        if not self.color:
            self.color = '0,0,100' if self.track_type == 'bigwig' else '0,100,0'
        
        # Set default visibility if not specified
        if not self.visibility:
            self.visibility = 'full' if self.track_type == 'bigwig' else 'pack'
        
        # Base track parameters (compatible with both browsers)
        track_params = {
            'name': track_name,
            'source': self.file_path,
            'visibility': self.visibility,
            'tracktype': self.track_type,
            'short_label': self.short_label,
            'long_label': self.long_label,
            'color': self.color,
            'priority': self.priority
        }
        
        # Add parent if provided
        if parent:
            track_params['parent'] = parent
        
        # Add track type specific defaults
        if self.track_type == 'bigwig':
            type_defaults = {
                'autoScale': 'on',
                'alwaysZero': 'on'
            }
            # Only add defaults if not overridden in additional_params
            for key, value in type_defaults.items():
                if key not in self.additional_params:
                    track_params[key] = value
        
        # Add Ensembl-specific parameters if requested
        if ensembl_compatible:
            # Add type and format parameters required by Ensembl
            if self.track_type == 'bigwig':
                track_params['type'] = 'bigWig'
                track_params['format'] = 'bigWig'
            elif self.track_type == 'bigbed':
                track_params['type'] = 'bigBed'
                track_params['format'] = 'bigBed'
        
        # Add any additional parameters
        for key, value in self.additional_params.items():
            track_params[key] = value
        
        return track_params


def parse_sample_sheet(file_path: str) -> Dict[str, SampleMetadata]:
    """
    Parse the sample sheet file.
    
    Args:
        file_path: Path to the sample sheet CSV file
        
    Returns:
        Dictionary mapping sample IDs to SampleMetadata objects
    """
    samples: Dict[str, SampleMetadata] = {}
    
    with open(file_path, 'r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        
        # Validate required columns
        required_columns = {'sample_id'}
        missing_columns = required_columns - set(reader.fieldnames)
        if missing_columns:
            raise ValueError(f"Sample sheet is missing required columns: {', '.join(missing_columns)}")
        
        for row in reader:
            # Skip empty rows
            if not row['sample_id'].strip():
                continue
                
            sample_id = row['sample_id']
            samples[sample_id] = SampleMetadata(row)
    
    return samples


def find_track_files(
    file_patterns: List[str],
    track_type: str,
    sample_pattern: Optional[Pattern] = None,
    annotation_pattern: Optional[Pattern] = None
) -> List[TrackFile]:
    """
    Find track files matching the provided glob patterns.
    
    Args:
        file_patterns: List of glob patterns to match files
        track_type: Type of track (bigwig or bigbed)
        sample_pattern: Regex pattern to extract sample IDs
        annotation_pattern: Regex pattern to extract annotation types
        
    Returns:
        List of TrackFile objects
    """
    track_files: List[TrackFile] = []
    
    for pattern in file_patterns:
        expanded_pattern = os.path.expanduser(pattern)
        files = glob.glob(expanded_pattern)
        
        if not files:
            logger.warning(f"No {track_type} files found matching pattern: {pattern}")
            continue
        
        logger.info(f"Found {len(files)} {track_type} files matching pattern: {pattern}")
        
        for file_path in files:
            track_file = TrackFile(
                file_path=file_path,
                track_type=track_type,
                sample_pattern=sample_pattern,
                annotation_pattern=annotation_pattern
            )
            track_files.append(track_file)
    
    return track_files


def apply_sample_metadata(
    track_files: List[TrackFile],
    sample_metadata: Dict[str, SampleMetadata]
) -> None:
    """
    Apply sample metadata to track files.
    
    Args:
        track_files: List of TrackFile objects
        sample_metadata: Dictionary mapping sample IDs to SampleMetadata objects
    """
    for track_file in track_files:
        if track_file.sample_id in sample_metadata:
            track_file.apply_metadata(sample_metadata[track_file.sample_id])
        else:
            logger.warning(f"No metadata found for sample ID: {track_file.sample_id} ({track_file.filename})")


def group_track_files_by_sample(
    track_files: List[TrackFile]
) -> Dict[str, List[TrackFile]]:
    """
    Group track files by sample ID.
    
    Args:
        track_files: List of TrackFile objects
        
    Returns:
        Dictionary mapping sample IDs to lists of TrackFile objects
    """
    grouped: Dict[str, List[TrackFile]] = {}
    
    for track_file in track_files:
        sample_id = track_file.sample_id
        if sample_id not in grouped:
            grouped[sample_id] = []
        grouped[sample_id].append(track_file)
    
    return grouped


def group_bigbed_files_by_annotation(
    track_files: List[TrackFile]
) -> Dict[str, List[TrackFile]]:
    """
    Group BigBed track files by annotation type.
    
    Args:
        track_files: List of TrackFile objects
        
    Returns:
        Dictionary mapping annotation types to lists of TrackFile objects
    """
    grouped: Dict[str, List[TrackFile]] = {}
    
    for track_file in track_files:
        annotation_type = track_file.annotation_type or 'default'
        if annotation_type not in grouped:
            grouped[annotation_type] = []
        grouped[annotation_type].append(track_file)
    
    return grouped


def convert_chromosome_names(
    file_path: str, 
    output_dir: Path,
    conversion_type: str = 'ucsc_to_ensembl'
) -> str:
    """
    Convert chromosome names in BigWig or BigBed files.
    
    Args:
        file_path: Path to the track file
        output_dir: Directory to write the converted file
        conversion_type: Type of conversion ('ucsc_to_ensembl' or 'ensembl_to_ucsc')
        
    Returns:
        Path to the converted file
    """
    # This is a placeholder for a real conversion function
    # In a real implementation, you would use tools like UCSC's wigToBigWig or bedToBigBed
    # with chromosome name mapping to convert between naming conventions
    
    # For demonstration purposes, we'll just copy the file and log a message
    basename = os.path.basename(file_path)
    output_path = output_dir / f"converted_{basename}"
    
    # In a real implementation, this would be a conversion operation
    logger.info(f"Would convert chromosome names in {file_path} from {conversion_type}")
    
    # Return the (imaginary) converted file path
    return str(output_path)


def create_data_hub(
    track_files: List[TrackFile],
    hub_name: str,
    genome: str,
    output_dir: Path,
    hub_email: Optional[str] = None,
    hub_url: Optional[str] = None,
    hub_description: Optional[str] = None,
    ensembl_compatible: bool = False,
    convert_chrom_names: bool = False
) -> None:
    """
    Create a track hub for data files.
    
    Args:
        track_files: List of TrackFile objects (all should be of the same track_type)
        hub_name: Name for the hub
        genome: Genome assembly (e.g., 'hg38')
        output_dir: Directory to write the hub files
        hub_email: Optional contact email for the hub
        hub_url: Optional URL where the hub will be hosted
        hub_description: Optional description for the hub
        ensembl_compatible: If True, include parameters for Ensembl compatibility
        convert_chrom_names: If True, convert chromosome names for Ensembl compatibility
    """
    if not track_files:
        logger.warning(f"No track files provided for hub {hub_name}")
        return
        
    hub_dir = output_dir / hub_name
    hub_dir.mkdir(exist_ok=True, parents=True)
    
    # Determine track type (all files should be the same type)
    track_type = track_files[0].track_type
    
    # Initialize hub
    hub = trackhub.Hub(
        hub_name,
        short_label=hub_name,
        long_label=f"{hub_name} - {track_type.upper()} Tracks",
        email=hub_email,
        descriptionUrl=hub_description
    )
    
    # Initialize genome
    # Ensure genome name is correctly formatted for the browser
    if ensembl_compatible and genome.startswith('hg'):
        # Convert UCSC genome names to Ensembl equivalents if needed
        # hg38 -> GRCh38, hg19 -> GRCh37
        ensembl_genome = genome.replace('hg38', 'GRCh38').replace('hg19', 'GRCh37')
        logger.info(f"Converting genome name from {genome} to {ensembl_genome} for Ensembl compatibility")
        genome_obj = trackhub.Genome(ensembl_genome)
    else:
        genome_obj = trackhub.Genome(genome)
    
    hub.add_genome(genome_obj)
    
    # Create a trackdb
    trackdb = trackhub.TrackDb()
    genome_obj.add_trackdb(trackdb)
    
    # Group files by sample
    samples = group_track_files_by_sample(track_files)
    
    # Create a composite track for each sample
    for sample_id, sample_files in samples.items():
        if not sample_files:
            continue
        
        # For Ensembl compatibility, we may need individual tracks instead of composites
        if ensembl_compatible:
            # Ensembl has limited support for composite tracks, so add tracks directly to trackdb
            for track_file in sample_files:
                track_params = track_file.get_track_params(ensembl_compatible=True)
                track = trackhub.Track(**track_params)
                trackdb.add_tracks(track)
        else:
            # UCSC-style composite tracks
            composite = trackhub.CompositeTrack(
                name=f"sample_{sample_id}",
                short_label=f"{sample_id}",
                long_label=f"Tracks for {sample_id}",
                visibility="full"
            )
            trackdb.add_tracks(composite)
            
            # Add tracks for each file in the sample
            for track_file in sample_files:
                track_params = track_file.get_track_params(parent=composite)
                track = trackhub.Track(**track_params)
                composite.add_tracks(track)
    
    # Write the hub to disk
    trackhub.upload.write_trackdb(hub, str(hub_dir), url=hub_url)
    logger.info(f"{track_type.upper()} hub written to {hub_dir}")
    
    # Generate URL for the hub
    if hub_url:
        hub_txt_url = f"{hub_url.rstrip('/')}/{hub_name}/hub.txt"
        logger.info(f"Hub URL: {hub_txt_url}")


def create_annotation_hub(
    track_files: List[TrackFile],
    annotation_type: str,
    hub_name: str,
    genome: str,
    output_dir: Path,
    hub_email: Optional[str] = None,
    hub_url: Optional[str] = None,
    hub_description: Optional[str] = None,
    ensembl_compatible: bool = False,
    convert_chrom_names: bool = False
) -> None:
    """
    Create a track hub for annotation files of a specific type.
    
    Args:
        track_files: List of TrackFile objects for this annotation type
        annotation_type: The annotation type/tool/parameter set
        hub_name: Base name for the hub
        genome: Genome assembly (e.g., 'hg38')
        output_dir: Directory to write the hub files
        hub_email: Optional contact email for the hub
        hub_url: Optional URL where the hub will be hosted
        hub_description: Optional description for the hub
        ensembl_compatible: If True, include parameters for Ensembl compatibility
        convert_chrom_names: If True, convert chromosome names for Ensembl compatibility
    """
    if not track_files:
        logger.warning(f"No track files provided for annotation type {annotation_type}")
        return
        
    # Create a specific hub name for this annotation type
    specific_hub_name = f"{hub_name}_{annotation_type}"
    hub_dir = output_dir / specific_hub_name
    hub_dir.mkdir(exist_ok=True, parents=True)
    
    # Initialize hub
    hub = trackhub.Hub(
        specific_hub_name,
        short_label=f"{annotation_type}",
        long_label=f"{hub_name} - {annotation_type} Annotations",
        email=hub_email,
        descriptionUrl=hub_description
    )
    
    # Initialize genome with browser-specific name if needed
    if ensembl_compatible and genome.startswith('hg'):
        ensembl_genome = genome.replace('hg38', 'GRCh38').replace('hg19', 'GRCh37')
        logger.info(f"Converting genome name from {genome} to {ensembl_genome} for Ensembl compatibility")
        genome_obj = trackhub.Genome(ensembl_genome)
    else:
        genome_obj = trackhub.Genome(genome)
    
    hub.add_genome(genome_obj)
    
    # Create a trackdb
    trackdb = trackhub.TrackDb()
    genome_obj.add_trackdb(trackdb)
    
    # Add a README track with information about this annotation type
    # Note: Ensembl may not support HTML tracks as UCSC does
    if not ensembl_compatible:
        readme_text = (
            f"# {annotation_type} Annotations\n\n"
            f"This track hub contains annotations generated using {annotation_type}.\n"
            f"These tracks should be viewed alongside the corresponding data tracks.\n"
        )
        
        readme_track = trackhub.HTMLTrack(
            name="readme",
            short_label="README",
            long_label=f"README for {annotation_type} Annotations",
            html=readme_text,
            visibility="hide"
        )
        trackdb.add_tracks(readme_track)
    
    # Group files by sample
    samples = group_track_files_by_sample(track_files)
    
    # Create tracks based on browser compatibility
    for sample_id, sample_files in samples.items():
        if not sample_files:
            continue
        
        # For Ensembl compatibility, we may need individual tracks instead of composites
        if ensembl_compatible:
            # Ensembl has limited support for composite tracks, so add tracks directly to trackdb
            for track_file in sample_files:
                track_params = track_file.get_track_params(ensembl_compatible=True)
                track = trackhub.Track(**track_params)
                trackdb.add_tracks(track)
        else:
            # UCSC-style composite tracks
            composite = trackhub.CompositeTrack(
                name=f"sample_{sample_id}",
                short_label=f"{sample_id}",
                long_label=f"{annotation_type} annotations for {sample_id}",
                visibility="pack"
            )
            trackdb.add_tracks(composite)
            
            # Add tracks for each file in the sample
            for track_file in sample_files:
                track_params = track_file.get_track_params(parent=composite)
                track = trackhub.Track(**track_params)
                composite.add_tracks(track)
    
    # Write the hub to disk
    trackhub.upload.write_trackdb(hub, str(hub_dir), url=hub_url)
    logger.info(f"Annotation hub for {annotation_type} written to {hub_dir}")
    
    # Generate URL for the hub
    if hub_url:
        hub_txt_url = f"{hub_url.rstrip('/')}/{specific_hub_name}/hub.txt"
        logger.info(f"Hub URL: {hub_txt_url}")


def create_sample_sheet_template(output_file: str) -> None:
    """
    Create a template sample sheet with example rows.
    
    Args:
        output_file: Path to write the template sample sheet to
    """
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = [
            'sample_id', 'short_label', 'long_label', 
            'color', 'visibility', 'priority', 'additional_params'
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        # Example rows
        writer.writerow({
            'sample_id': 'Sample1',
            'short_label': 'Sample 1',
            'long_label': 'Sample 1 - Wild Type',
            'color': '0,0,100',
            'visibility': 'full',
            'priority': '1',
            'additional_params': '{"autoScale":"on","alwaysZero":"on"}'
        })
        writer.writerow({
            'sample_id': 'Sample2',
            'short_label': 'Sample 2',
            'long_label': 'Sample 2 - Treatment',
            'color': '100,0,0',
            'visibility': 'full',
            'priority': '2',
            'additional_params': '{"autoScale":"on","alwaysZero":"on"}'
        })
    
    logger.info(f"Sample sheet template written to {output_file}")


def main() -> None:
    """
    Main function to parse arguments and create track hubs.
    """
    parser = argparse.ArgumentParser(
        description="Create UCSC and Ensembl compatible track hubs from globbed files with sample metadata",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Command mode subparsers
    subparsers = parser.add_subparsers(dest='command', help='Command to execute')
    
    # Create template sample sheet
    template_parser = subparsers.add_parser(
        'template', 
        help='Create a template sample sheet'
    )
    template_parser.add_argument(
        'output_file',
        help='Output file for the template sample sheet'
    )
    
    # Create hubs from globbed files with sample metadata
    create_parser = subparsers.add_parser(
        'create',
        help='Create track hubs from globbed files with sample metadata'
    )
    create_parser.add_argument(
        '--sample-sheet',
        help='Path to the optional sample sheet CSV file with metadata'
    )
    create_parser.add_argument(
        '--bigwig',
        nargs='+',
        help='Glob patterns for BigWig files (e.g., "data/*.bw")'
    )
    create_parser.add_argument(
        '--bigbed',
        nargs='+',
        help='Glob patterns for BigBed files (e.g., "annotations/*.bb")'
    )
    create_parser.add_argument(
        '--hub-name',
        required=True,
        help='Base name for the hubs (e.g., "RiboSeq")'
    )
    create_parser.add_argument(
        '--genome',
        required=True,
        help='Genome assembly (e.g., "hg38")'
    )
    create_parser.add_argument(
        '--output-dir',
        required=True,
        help='Directory to write the hubs to'
    )
    create_parser.add_argument(
        '--sample-regex',
        help='Regex pattern with named group "sample_id" to extract sample IDs from filenames'
    )
    create_parser.add_argument(
        '--annotation-regex',
        help='Regex pattern with named group "annotation_type" to identify annotation types'
    )
    create_parser.add_argument(
        '--email',
        help='Contact email for the hubs'
    )
    create_parser.add_argument(
        '--hub-url',
        help='URL where the hubs will be hosted'
    )
    create_parser.add_argument(
        '--hub-description',
        help='Path to a file containing hub description in HTML format'
    )
    create_parser.add_argument(
        '--ensembl-compatible',
        action='store_true',
        help='Make track hubs compatible with Ensembl'
    )
    create_parser.add_argument(
        '--convert-chrom-names',
        action='store_true',
        help='Convert chromosome names between UCSC and Ensembl formats (requires liftOver tools)'
    )
    
    args = parser.parse_args()
    
    # Handle template command
    if args.command == 'template':
        create_sample_sheet_template(args.output_file)
        return
        
    # Handle create command
    if args.command == 'create':
        # Check that at least one of bigwig or bigbed is specified
        if not args.bigwig and not args.bigbed:
            logger.error("At least one of --bigwig or --bigbed must be specified")
            return
        
        # Parse hub description if provided
        hub_description = None
        if args.hub_description and os.path.exists(args.hub_description):
            with open(args.hub_description, 'r') as f:
                hub_description = f.read()
        
        # Create output directory
        output_dir = Path(args.output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        
        # Compile regex patterns if provided
        sample_pattern = None
        if args.sample_regex:
            sample_pattern = re.compile(args.sample_regex)
        
        annotation_pattern = None
        if args.annotation_regex:
            annotation_pattern = re.compile(args.annotation_regex)
        
        # Load sample metadata if provided
        sample_metadata: Dict[str, SampleMetadata] = {}
        if args.sample_sheet:
            sample_metadata = parse_sample_sheet(args.sample_sheet)
            logger.info(f"Loaded metadata for {len(sample_metadata)} samples from {args.sample_sheet}")
        
        # Find and process BigWig files if specified
        if args.bigwig:
            bigwig_files = find_track_files(
                args.bigwig, 'bigwig', sample_pattern, None
            )
            
            if bigwig_files:
                # Apply sample metadata if available
                if sample_metadata:
                    apply_sample_metadata(bigwig_files, sample_metadata)
                
                # Create data hub
                create_data_hub(
                    track_files=bigwig_files,
                    hub_name=f"{args.hub_name}_Data",
                    genome=args.genome,
                    output_dir=output_dir,
                    hub_email=args.email,
                    hub_url=args.hub_url,
                    hub_description=hub_description,
                    ensembl_compatible=args.ensembl_compatible,
                    convert_chrom_names=args.convert_chrom_names
                )
       
       # Find and process BigBed files if specified
    if args.bigbed:
        bigbed_files = find_track_files(
            args.bigbed, 'bigbed', sample_pattern, annotation_pattern
        )
        
        if bigbed_files:
            # Apply sample metadata if available
            if sample_metadata:
                apply_sample_metadata(bigbed_files, sample_metadata)
            
            # Group by annotation type
            annotation_groups = group_bigbed_files_by_annotation(bigbed_files)
            
            # Create a hub for each annotation type
            for annotation_type, annotation_files in annotation_groups.items():
                create_annotation_hub(
                    track_files=annotation_files,
                    annotation_type=annotation_type,
                    hub_name=f"{args.hub_name}_Annotation",
                    genome=args.genome,
                    output_dir=output_dir,
                    hub_email=args.email,
                    hub_url=args.hub_url,
                    hub_description=hub_description,
                    ensembl_compatible=args.ensembl_compatible,
                    convert_chrom_names=args.convert_chrom_names
                )
            
            logger.info(f"Created {len(annotation_groups)} annotation hubs")
       
       # Provide instructions for using the hubs
    print("\nTrack Hub Creation Complete!")
       print(f"Track hubs have been created in: {output_dir}")
       
       if args.hub_url:
           if args.bigwig:
               data_hub_url = f"{args.hub_url.rstrip('/')}/{args.hub_name}_Data/hub.txt"
               print(f"\nTo add the data hub to UCSC Genome Browser:")
               print(f"1. Go to https://genome.ucsc.edu/cgi-bin/hgHubConnect")
               print(f"2. Click 'My Hubs' tab")
               print(f"3. Paste this URL: {data_hub_url}")
               print(f"4. Click 'Add Hub'")
               
               if args.ensembl_compatible:
                   print(f"\nTo add the data hub to Ensembl:")
                   print(f"1. Go to https://www.ensembl.org")
                   print(f"2. Navigate to a species page")
                   print(f"3. Click 'Add your data' (in the left menu)")
                   print(f"4. Select 'Add Track Hub'")
                   print(f"5. Paste this URL: {data_hub_url}")
           
           if args.bigbed and annotation_groups:
               print("\nAnnotation hubs:")
               for annotation_type in annotation_groups.keys():
                   annotation_hub_url = f"{args.hub_url.rstrip('/')}/{args.hub_name}_Annotation_{annotation_type}/hub.txt"
                   print(f"- {annotation_type}: {annotation_hub_url}")
       else:
           print("\nTo use these track hubs, you need to:")
           print("1. Host the hub directory on a web server")
           print("2. Use the URL to the hub.txt file when adding the hub to UCSC Genome Browser or Ensembl")
       
       return
   
   # If no command specified, show help
   parser.print_help()


if __name__ == "__main__":
   main()