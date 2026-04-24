"""This module contains DTO for project configuration."""
import enum
import yaml
import typing as tp
from pydantic import BaseModel


class DatasetTag(str, enum.Enum):
    # Acquisition
    synthetic = "synthetic"
    nontryptic = "nontryptic"
    timstof = "timstof"
    waters = "waters"
    sciex = "sciex"
    agilent = "agilent"
    astral = "astral"
    # Labeling strategy
    silac = "silac"
    tmt = "tmt"
    # Common modifications
    carbamidomethylation = "carbamidomethylation"
    deamidation = "deamidation"
    oxidation = "oxidation"
    acetylation = "acetylation"
    methylation = "methylation"
    dimethylation = "dimethylation"
    trimethylation = "trimethylation"
    phosphorylation = "phosphorylation"
    ammonia_loss = "ammonia_loss"
    # Other modifications
    sodium_adduct = "sodium_adduct"
    formylation = "formylation"
    formaldehyde = "formaldehyde"
    acetaldehyde = "acetaldehyde"
    propionyl = "propionyl"
    butyryl = "butyryl"
    crotonylation = "crotonylation"
    hydroxyisobutyryl = "hydroxyisobutyryl"
    biotinylation = "biotinylation"
    ubiquitination = "ubiquitination"
    succinyl = "succinyl"
    malonylation = "malonylation"
    glutarylation = "glutarylation"
    citrullination = "citrullination"
    nitro = "nitro"
    diethyl = "diethyl"
    nickel_2h_substitution = "nickel_2h_substitution" # # Replacement of 2 protons by nickel
    iron_2h_substitution = "iron_2h_substitution" # Replacement of 2 protons by iron
    iron_3h_substitution = "iron_3h_substitution" # Replacement of 3 protons by iron


class Instrument(str, enum.Enum):
    # MS instrument types distinguished by MSGF+
    LCQ = "LCQ" # Linear Quadrupole Ion Trap
    LTQ = "LTQ" # Linear Trap Quadrupole
    Orbitrap = "Orbitrap" # Orbitrap Elite, Fusion, Lumos
    FTICR = "FTICR" # Fourier Transform Ion Cyclotron Resonance
    Lumos = "Lumos" # a high-end version of the Orbitrap Fusion
    TOF = "TOF" # Bruker, Sciex TripleTOF, Agilent TOF
    QExactive = "QExactive" # Hybrid Quadrupole-Orbitrap (QE, QE Plus / HF)


class FragmentationMethod(str, enum.Enum):
    not_given = "not_given" # TODO: mb rename! "as written in the spectrum or CID if no info"
    CID = "CID"
    ETD = "ETD"
    HCD = "HCD"
    UVPD = "UVPD"


class DataDownloadConfig(BaseModel):
    # Dataset id in MS repository 
    dset_id: str
    # Download links (alternative to dataset id)
    links: tp.Optional[tp.List[str]] = None
    # Extension of raw spectra files to download (lowercase)
    ext: str = ".raw"
    # Number of raw files to download.
    # If null, all matching files are downloaded
    n_files: tp.Optional[int] = None
    # Keywords to select files to download.
    # Only filenames containing all the keywords are downloaded
    keywords: tp.Optional[tp.List[str]] = []
    # Keywords to exclude files from downloading.
    exclude_keywords: tp.Optional[tp.List[str]] = []


class DBSearchConfig(BaseModel):
    # Path to the protein database file in FASTA format.
    database_path: tp.Optional[tp.Union[str, tp.List[str]]] = None
    # For synthetic peptides: path to dir with corresponding pool proteomes.
    # Will be appened to PROTEOMES_DIR.
    pool_proteomes_dir: tp.Optional[str] = None
    # Extension of spectra files used for the search
    ext: str = ".mzml"
    # Number of splits (for running DB search with large DB).
    n_db_splits: int = 1
    # Additional params for MSGF+
    instrument: Instrument = Instrument.Orbitrap # TODO: what is the default?
    fragmentation: FragmentationMethod = FragmentationMethod.HCD # TODO: set "CID" or "-" as default?
    # Additional (optional) search params. Only needs to be passed 
    # for non-default params. Can be an emtpy dict. 
    search_params: tp.Dict[str, tp.Any]


class RescoringConfig(BaseModel):
    # Additional feature prediction params for MSBooster.
    # Can be an empty dict.
    feat_pred_params: tp.Dict[str, tp.Any]
    # Max q-value threshold for rescored DB search results
    q_val_threshold: tp.Optional[float] = 0.01
    # If specified, filter only peptides from the pools in the list.
    # Only relevant for searches with a single protein database file.
    # For searches with pools from pool_proteomes_dir, always filter by the corresponding pool.
    filter_pools: tp.Optional[tp.List[str]] = None


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
