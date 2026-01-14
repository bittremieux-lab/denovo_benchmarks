#!/usr/bin/env python3
"""
Migration script to convert outputs from OLD to NEW versioned structure.

OLD structure:
  outputs/{dataset}/{algorithm}_output.csv

NEW structure:
  outputs/{algorithm}/{version}/{dataset}_output.csv
  outputs/{algorithm}/latest/ -> {version}  (symlink)

This script COPIES files (does not move/delete) for safety.
Use --cleanup flag to remove old structure after successful migration.

Usage:
  python migrate_to_versioned_structure.py --outputs-dir outputs/
  python migrate_to_versioned_structure.py --outputs-dir tests/mock_outputs/ --dry-run
  python migrate_to_versioned_structure.py --outputs-dir outputs/ --cleanup
"""

import argparse
import shutil
from pathlib import Path
from typing import Dict, List, Tuple

from version_utils import get_all_algorithms_versions, create_latest_symlink


def discover_old_structure_files(outputs_dir: Path) -> List[Tuple[Path, str, str]]:
    """
    Discover all output files in the OLD structure.
    
    Args:
        outputs_dir: Base outputs directory
        
    Returns:
        List of tuples: (file_path, algorithm_name, dataset_name)
    """
    files = []
    
    if not outputs_dir.exists():
        print(f"Warning: Directory {outputs_dir} does not exist")
        return files
    
    for dataset_dir in outputs_dir.iterdir():
        if not dataset_dir.is_dir():
            continue
        
        dataset_name = dataset_dir.name
        
        # Find all output files in this dataset
        for output_file in dataset_dir.glob("*_output*.csv"):
            filename = output_file.stem
            
            # Parse algorithm name from filename
            if filename.endswith('_output_augmented'):
                algorithm = filename.replace('_output_augmented', '')
                suffix = '_output_augmented'
            elif filename.endswith('_output'):
                algorithm = filename.replace('_output', '')
                suffix = '_output'
            else:
                print(f"Warning: Unexpected file format: {output_file.name}")
                continue
            
            files.append((output_file, algorithm, dataset_name, suffix))
    
    return files


def migrate_outputs(
    outputs_dir: Path,
    algorithms_dir: Path = Path("algorithms"),
    dry_run: bool = False,
    cleanup: bool = False
) -> None:
    """
    Migrate outputs from OLD to NEW structure by copying files.
    
    Args:
        outputs_dir: Base outputs directory to migrate
        algorithms_dir: Directory containing algorithm versions.log files
        dry_run: If True, only show what would be done without making changes
        cleanup: If True, remove old structure after successful migration
    """
    print("="*70)
    print("MIGRATION: OLD STRUCTURE -> NEW VERSIONED STRUCTURE")
    print("="*70)
    
    if dry_run:
        print("DRY RUN MODE - No changes will be made\n")
    
    # Get all algorithm versions
    print("1. Reading algorithm versions from versions.log files...")
    versions_map = get_all_algorithms_versions(algorithms_dir)
    
    if not versions_map:
        print("ERROR: No algorithms with versions.log found!")
        print(f"   Check {algorithms_dir} directory")
        return
    
    for algo, version in sorted(versions_map.items()):
        print(f"   {algo:20s} -> {version}")
    
    # Discover files in old structure
    print(f"\n2. Discovering files in OLD structure: {outputs_dir}")
    old_files = discover_old_structure_files(outputs_dir)
    
    if not old_files:
        print("   WARNING: No files found in OLD structure. Nothing to migrate.")
        return
    
    print(f"   Found {len(old_files)} files to migrate")
    
    # Group by algorithm
    files_by_algo: Dict[str, List] = {}
    for file_path, algorithm, dataset, suffix in old_files:
        if algorithm not in files_by_algo:
            files_by_algo[algorithm] = []
        files_by_algo[algorithm].append((file_path, dataset, suffix))
    
    # Migrate each algorithm
    print("\n3. Copying files to NEW structure...")
    
    copied_files = 0
    created_dirs = set()
    created_symlinks = set()
    
    for algorithm in sorted(files_by_algo.keys()):
        print(f"\n   Algorithm: {algorithm}")
        
        # Get version for this algorithm
        if algorithm not in versions_map:
            print("      WARNING: No version found in versions.log, skipping")
            continue
        
        version = versions_map[algorithm]
        print(f"      Version: {version}")
        
        # Create new directory structure
        new_algo_dir = outputs_dir / algorithm / version
        
        if not dry_run:
            new_algo_dir.mkdir(parents=True, exist_ok=True)
        
        created_dirs.add(new_algo_dir)
        
        # Copy files
        for file_path, dataset, suffix in files_by_algo[algorithm]:
            # New filename: {dataset}_output.csv or {dataset}_output_augmented.csv
            new_filename = f"{dataset}{suffix}.csv"
            new_file_path = new_algo_dir / new_filename
            
            if dry_run:
                print(f"      [would copy] {file_path.name} -> {new_file_path.relative_to(outputs_dir)}")
            else:
                shutil.copy2(str(file_path), str(new_file_path))
                print(f"      Copied: {file_path.name} -> {new_file_path.relative_to(outputs_dir)}")
            
            copied_files += 1
        
        # Create symlink
        if not dry_run:
            try:
                create_latest_symlink(outputs_dir, algorithm, version)
                created_symlinks.add(algorithm)
                print(f"      Symlink: latest -> {version}")
            except Exception as e:
                print(f"      WARNING: Failed to create symlink: {e}")
        else:
            print(f"      [would create] latest -> {version}")
            created_symlinks.add(algorithm)
    
    # Clean up old structure (if requested)
    print("\n4. Cleaning up old structure...")
    
    removed_files = 0
    removed_dirs = 0
    
    if cleanup and not dry_run:
        # Remove the old files and empty directories
        print("   Cleanup mode - Removing old structure")
        
        # Remove old files
        for file_path, _, _, _ in old_files:
            if file_path.exists():
                file_path.unlink()
                print(f"   Removed: {file_path}")
                removed_files += 1
        
        # Remove empty directories
        for dataset_dir in outputs_dir.iterdir():
            if dataset_dir.is_dir() and not any(dataset_dir.iterdir()):
                dataset_dir.rmdir()
                print(f"   Removed empty directory: {dataset_dir.name}")
                removed_dirs += 1
    elif cleanup and dry_run:
        # Show what would be cleaned up
        print("   [would remove] old files and directories")
        
        old_dataset_dirs = set()
        for file_path, _, _, _ in old_files:
            old_dataset_dirs.add(file_path.parent)
            removed_files += 1
        removed_dirs = len(old_dataset_dirs)
    else:
        print("   Old structure preserved (use --cleanup to remove)")
    
    # Summary
    print("\n" + "="*70)
    print("MIGRATION SUMMARY")
    print("="*70)
    
    if dry_run:
        print(f"   Would copy:    {copied_files} files")
        print(f"   Would create:  {len(created_dirs)} version directories")
        print(f"   Would create:  {len(created_symlinks)} symlinks")
        if cleanup:
            print(f"   Would remove:  {removed_files} old files")
            print(f"   Would remove:  {removed_dirs} empty directories")
        print("\n   Run without --dry-run to perform migration")
    else:
        print(f"   Copied:      {copied_files} files")
        print(f"   Created:     {len(created_dirs)} version directories")
        print(f"   Created:     {len(created_symlinks)} symlinks")
        if cleanup:
            print(f"   Removed:     {removed_files} old files")
            print(f"   Removed:     {removed_dirs} empty directories")
        else:
            print("   Old files preserved (use --cleanup to remove)")
        print("\n   Migration complete!")
    
    # Show new structure
    print("\n" + "="*70)
    print("NEW STRUCTURE")
    print("="*70)
    
    for algorithm in sorted(files_by_algo.keys()):
        if algorithm in versions_map:
            version = versions_map[algorithm]
            print(f"   {outputs_dir.name}/{algorithm}/")
            print(f"      {version}/")
            
            # Count files
            if not dry_run:
                new_algo_dir = outputs_dir / algorithm / version
                if new_algo_dir.exists():
                    n_files = len(list(new_algo_dir.glob("*.csv")))
                    print(f"         ({n_files} CSV files)")
            
            if algorithm in created_symlinks:
                print(f"      latest/ -> {version}")


def main():
    parser = argparse.ArgumentParser(
        description="Migrate outputs from OLD to NEW versioned structure",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Dry run on real outputs
  python migrate_to_versioned_structure.py --outputs-dir outputs/ --dry-run
  
  # Test on mock data
  python migrate_to_versioned_structure.py --outputs-dir tests/mock_outputs/ --dry-run
  
  # Actually perform migration
  python migrate_to_versioned_structure.py --outputs-dir outputs/
        """
    )
    
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        default=Path("outputs"),
        help="Base outputs directory to migrate (default: outputs/)"
    )
    
    parser.add_argument(
        "--algorithms-dir",
        type=Path,
        default=Path("algorithms"),
        help="Directory containing algorithm versions.log files (default: algorithms/)"
    )
    
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without making changes"
    )
    
    parser.add_argument(
        "--cleanup",
        action="store_true",
        help="Remove old structure after successful migration (by default, old files are preserved)"
    )
    
    args = parser.parse_args()
    
    migrate_outputs(
        outputs_dir=args.outputs_dir,
        algorithms_dir=args.algorithms_dir,
        dry_run=args.dry_run,
        cleanup=args.cleanup
    )


if __name__ == "__main__":
    main()
