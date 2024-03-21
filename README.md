# Benchmarking de novo* peptide sequencing algorithms

## Requirements

```
pip install -r requirements.txt
```

## Usage

1. Place the dataset in the `datasets` folder (or use `datasets/sample_data`).

2. Execute the `run.sh` script from the command line, passing the path to your dataset as an argument:

```
./run.sh path/to/your/data
```

This script will create and run computations within Docker containers, collect results, and evaluate the algorithms located in the `algorithms` folder.

## Output

- Algorithm predictions are collected in the `outputs/` folder in separate files named `algorithm_outputs.csv`.
- Evaluation results are stored in `metrics.csv`.
