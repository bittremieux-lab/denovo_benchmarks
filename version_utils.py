"""
Version management utilities for the de novo benchmarking system.
"""

import yaml
from pathlib import Path
from typing import Dict, List, Optional


ALGORITHMS_DIR = Path("algorithms")


def parse_versions_log(algorithm_path: Path) -> List[Dict]:
    """
    Parse versions.log file from an algorithm directory.
    """
    versions_file = algorithm_path / "versions.log"
    
    if not versions_file.exists():
        raise FileNotFoundError(
            f"versions.log not found in {algorithm_path}. "
            f"Please create one based on algorithms/base/versions_template.log"
        )
    
    with open(versions_file, 'r') as f:
        versions = yaml.safe_load(f)
    
    if not versions or not isinstance(versions, list):
        raise ValueError(
            f"Invalid versions.log format in {algorithm_path}. "
            f"Expected a list of version entries."
        )
    
    return versions


def get_latest_version(algorithm_path: Path) -> str:
    """
    Get the most recent version from versions.log.
    """
    versions = parse_versions_log(algorithm_path)
    return versions[0]['container_version']


def get_all_algorithms_versions(algorithms_base_dir: Optional[Path] = None) -> Dict[str, str]:
    """
    Scan all algorithm directories and get their latest versions.
    """
    if algorithms_base_dir is None:
        algorithms_base_dir = ALGORITHMS_DIR
    
    versions_map = {}
    
    if not algorithms_base_dir.exists():
        return versions_map
    
    for algo_dir in algorithms_base_dir.iterdir():
        if not algo_dir.is_dir():
            continue
        
        # Skip base directory (contains templates)
        if algo_dir.name == "base":
            continue
        
        # Skip if no versions.log exists
        if not (algo_dir / "versions.log").exists():
            print(f"Warning: No versions.log found for {algo_dir.name}, skipping.")
            continue
        
        try:
            version = get_latest_version(algo_dir)
            versions_map[algo_dir.name] = version
        except Exception as e:
            print(f"Warning: Could not parse version for {algo_dir.name}: {e}")
            continue
    
    return versions_map


def create_latest_symlink(outputs_base_dir: Path, algorithm: str, version: str) -> None:
    """
    Create or update the 'latest' symlink for an algorithm.
    Creates: outputs/{algorithm}/latest -> {version}
    """
    algo_dir = outputs_base_dir / algorithm
    latest_link = algo_dir / "latest"
    version_dir = algo_dir / version
    
    # Ensure the version directory exists
    if not version_dir.exists():
        raise FileNotFoundError(
            f"Version directory {version_dir} does not exist. "
            f"Cannot create symlink to non-existent directory."
        )
    
    # Remove existing symlink if it exists
    if latest_link.is_symlink() or latest_link.exists():
        latest_link.unlink()
    
    # Create new symlink (relative path)
    latest_link.symlink_to(version, target_is_directory=True)
    print(f"Created symlink: {latest_link} -> {version}")
