#!/bin/bash

# Make pyenvs available
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"

# Run Casanovo to obtain initial starting sequences
pyenv activate casanovo
casanovo_output="outputs.mztab"

# Get dataset property tags
DSET_TAGS=$(python /algo/base/dataset_tags_parser.py --dataset "$@")
while IFS='=' read -r key value; do
    export "$key"="$value"
done <<< "$DSET_TAGS"

# Use tag variables to specify de novo algorithm
# for the particular dataset properties
if  [[ -v nontryptic && $nontryptic -eq 1 ]]; then
    # Run de novo algorithm on the input data
    echo "Using non-tryptic model."
    casanovo sequence -c casanovo_config.yml -o $casanovo_output "$@"/*.mgf --model ./casanovo_nontryptic.ckpt
else
    # Run de novo algorithm on the input data
    echo "Using general model."
    casanovo sequence -c casanovo_config.yml -o $casanovo_output "$@"/*.mgf
fi

# Change to spectralis pyenv
pyenv deactivate
pyenv activate spectralis

# Write the initial starting sequences to the input MGFs
spectralis_data_dir="./seq_data"
mkdir -p $spectralis_data_dir
python /algo/intermediate_mapper.py --mztab_path $casanovo_output --mgf_in_dir "$@" --mgf_out_dir $spectralis_data_dir

# Run Spectralis on every MGF file in the input directory
input_mgf_paths=$spectralis_data_dir/*.mgf
mkdir -p ./spectralis_outputs
spectralis_outputs=./spectralis_outputs/
for mgf_file in $input_mgf_paths; do
    file_name=$(basename -- "$mgf_file")
    # Something wrong with readgin spectra:
#     Traceback (most recent call last):
#   File "/opt/custom/.pyenv/versions/spectralis/bin/spectralis", line 8, in <module>
#     sys.exit(main())
#   File "/opt/custom/.pyenv/versions/3.7.17/envs/spectralis/lib/python3.7/site-packages/click/core.py", line 1157, in __call__
#     return self.main(*args, **kwargs)
#   File "/opt/custom/.pyenv/versions/3.7.17/envs/spectralis/lib/python3.7/site-packages/click/core.py", line 1078, in main
#     rv = self.invoke(ctx)
#   File "/opt/custom/.pyenv/versions/3.7.17/envs/spectralis/lib/python3.7/site-packages/click/core.py", line 1434, in invoke
#     return ctx.invoke(self.callback, **ctx.params)
#   File "/opt/custom/.pyenv/versions/3.7.17/envs/spectralis/lib/python3.7/site-packages/click/core.py", line 783, in invoke
#     return __callback(*args, **kwargs)
#   File "/opt/custom/.pyenv/versions/3.7.17/envs/spectralis/lib/python3.7/site-packages/spectralis/main.py", line 81, in main
#     spectralis.rescoring_from_mgf(mgf_path=input_path, out_path=output_path)
#   File "/opt/custom/.pyenv/versions/3.7.17/envs/spectralis/lib/python3.7/site-packages/spectralis/spectralis_master.py", line 510, in rescoring_from_mgf
#     _out = self._process_mgf(mgf_path)
#   File "/opt/custom/.pyenv/versions/3.7.17/envs/spectralis/lib/python3.7/site-packages/spectralis/spectralis_master.py", line 234, in _process_mgf
#     for spectrum in tqdm.tqdm(reader):
#   File "/opt/custom/.pyenv/versions/3.7.17/envs/spectralis/lib/python3.7/site-packages/tqdm/std.py", line 1181, in __iter__
#     for obj in iterable:
#   File "/opt/custom/.pyenv/versions/3.7.17/envs/spectralis/lib/python3.7/site-packages/pyteomics/auxiliary/file_helpers.py", line 178, in __next__
    spectralis --mode="rescoring" --config="spectralis_config.yml" --input_path="$mgf_file" --output_path="${spectralis_outputs}${file_name%.mgf}_spectralis_rescoring.csv"
done