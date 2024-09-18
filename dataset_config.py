"""This module contains DTO for project configuration."""
import enum
import yaml
import typing as tp
from pydantic import BaseModel


class DatasetTag(str, enum.Enum):
    synthetic = "synthetic"
    nontryptic = "nontryptic"
    timstof = "timstof"
    waters = "waters"
    sciex = "sciex"
    phosphorylation = "phosphorylation"


class DataDownloadConfig(BaseModel):
    # Dataset id in MS repository 
    dset_id: str
    # Download links (alternative to dataset id)
    links: tp.Optional[tp.List[str]] = None
    # Extension of raw spectra files to download (lowercase) # TODO: rename to raw_ext? 
    ext: str = ".raw"
    # Number of raw files to download.
    # If null, all matching files are downloaded
    n_files: tp.Optional[int] = None
    # Keywords to select files to download.
    # Only filenames containing all the keywords are downloaded
    keywords: tp.Optional[tp.List[str]] = []


class DBSearchConfig(BaseModel):
    # Path to the protein database file in FASTA format.
    database_path: tp.Union[str, tp.List[str]]
    # Extension of spectra files used for the search # TODO: rename to search_ext? 
    ext: str = ".mzml"
    # Number of splits (for running DB search with large DB).
    n_db_splits: int = 1
    # Additional (optional) search params. Only needs to be passed 
    # for non-default params. Can be an emtpy dict. 
    search_params: tp.Dict[str, tp.Any]


class RescoringConfig(BaseModel):
    # Additional feature prediction params for MSBooster.
    # Can be an empty dict.
    feat_pred_params: tp.Dict[str, tp.Any]
    # Max q-value threshold for rescored DB search results
    q_val_threshold: float = 0.01


class Config(BaseModel):
    """Keeps all dataset processing options."""
    # Dataset name
    name: str
    # Dataset properties tags
    tags: tp.List[DatasetTag] = []
    # Dataset description
    desc: str
    # Data downloading related properties
    download: DataDownloadConfig
    # Database search related properties
    db_search: DBSearchConfig
    # Rescoring related properties
    rescoring: RescoringConfig


def get_config(config_path: str) -> Config:
    """Parse .YAML file with project options and build options object.

    Parameters:
        config_path: Path to configuration .YAML file.

    Returns:
        Options serialized in object.
    """
    with open(config_path, "r") as yf:
        yml_file = yaml.safe_load(yf)
        config = Config.parse_obj(yml_file)
    return config
