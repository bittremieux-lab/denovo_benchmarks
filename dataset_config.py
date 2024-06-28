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


class DataDownloadConfig(BaseModel):
    # Dataset id in MS repository 
    dset_id: str
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
    # Additional (optional) search params. Only needs to be passed 
    # for non-default params. Can be an emtpy dict. 
    search_params: tp.Dict[str, tp.Any]


class RescoringConfig(BaseModel):
    # Additional feature prediction params for MSBooster.
    # Can be an empty dict.
    feat_pred_params: tp.Dict[str, tp.Any]
    # the format the spectra are provided in ("mzml", "RAW", "d")
    spectra_type: str # TODO: to enum?
    # search_results_type: str # "MSFragger"
    # the model used for fragment intensity prediction, e.g. "some model"
    intensity_model: str # TODO: to enum?
    # the model used for retention time prediction, e.g. "some model"
    irt_model: str # TODO: to enum?
    # TODO: Optional?
    massTolerance: float = 20
    unitMassTolerance: str = "ppm"
    # ce_range: tp.List[int] = [19, 50]
    # Max q-value threshold for rescored DB search results
    q_val_threshold: float = 0.01


# TODO: remove
# class SubModel(BaseModel):
#     foo: str = 'bar'
#     apple: int = 1


class Config(BaseModel):
    """Keeps all dataset processing options."""
    # Dataset name
    name: str
    # Dataset properties tags
    tags: tp.List[DatasetTag] = []
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
