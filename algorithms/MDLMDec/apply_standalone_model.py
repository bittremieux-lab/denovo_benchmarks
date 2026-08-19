from load_standalone_model import load_model
from loader_ import LoaderHF
import torch
from denovo_base.utils import Dict2dev
from tqdm import tqdm
import numpy as np
import pandas as pd
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

def fill_null_after_first_eos_token(intseq, EOS, NT):
    if len(intseq.shape) == 1:
        intseq = intseq[None]
    bs, sl = intseq.shape
    length = ((intseq == EOS)|(intseq == NT)).int().argmax(1)
    intseq[torch.arange(bs), length] = EOS
    mask = (length > 0)[:,None]
    index_array = torch.arange(sl)[None].repeat([bs, 1]).to(intseq.device)
    boolean_array = index_array > length[:, None]
    intseq[mask & boolean_array] = NT

    return intseq

def to_list_of_strings(intseq, reverse, revdic, EOS, NT):
    if len(intseq.shape) == 1:
        intseq = intseq[None]
    is_reverse = lambda x: x[::-1] if reverse else x
    return [
        is_reverse([
            revdic[int(n)] for n in m if n not in [NT, EOS]
        ])
        for m in intseq
    ]

def to_list_of_numbers(array, aaseqlist, reverse):
    if len(array.shape) == 1:
        array = array[None]
    is_reverse = lambda x: x[::-1] if reverse else x
    return [
        is_reverse([
            array[i, j].item() for j, n in enumerate(m)
        ])
        for i, m in enumerate(aaseqlist)
    ]

def main():
    model = load_model("./experiment_directory", device=device)

    loader = LoaderHF(
        "parquet",
        synonyms=[['I','L']],
        reverse=model.reverse,
        batch_size=100,
        num_workers=2,
    )

    pbar = tqdm(loader.dataloader['test'], total=loader.size//100)
    output = {
        'spectrum_id': [],
        'filename': [],
        'index': [],
        'prediction': [],
        'aa_prob': [],
        'pep_prob': [],
        'aa_prob_min': [],
        'pep_prob_min': [],
    }
    for batch in pbar:
        output['spectrum_id'].extend(batch['id'].tolist())
        output['filename'].extend(batch['file'].tolist())
        output['index'].extend(batch['index'].tolist())
        batch = Dict2dev(batch, device)
        outdict = model.predict_sequence(batch, save_p=True, save_x=True)
        # ['prediction', 'logits', 'x_save', 'p_save', 'aa_prob_min', 'pep_prob_min', 'aa_entropy', 'pep_entropy']
        prediction = outdict.pop('prediction')
        prediction = fill_null_after_first_eos_token(prediction, model.decoder.EOS, model.decoder.NT)
        pred_strings = to_list_of_strings(prediction, model.reverse, model.decoder.rev_outdict, model.decoder.EOS, model.decoder.NT)
        logits = outdict.pop('logits').softmax(-1).gather(-1, prediction[...,None]).squeeze()
        aaprob_list = to_list_of_numbers(logits, pred_strings, model.reverse)
        pepprob = [np.prod(lst).item() for lst in aaprob_list]
        aamin = outdict.pop('aa_prob_min')
        aamin_list = to_list_of_numbers(aamin, pred_strings, model.reverse)
        pepmin = outdict.pop("pep_prob_min")

        output['prediction'].extend(pred_strings)
        output['aa_prob'].extend(aaprob_list)
        output['aa_prob_min'].extend(aamin_list)
        output['pep_prob'].extend(pepprob)
        output['pep_prob_min'].extend(pepmin.cpu().tolist())

    pd.DataFrame(output).to_parquet("output.parquet")

if __name__ == '__main__':
    main()
