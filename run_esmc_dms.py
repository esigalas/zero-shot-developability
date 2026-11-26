# run_esmc_dms.py
import argparse
import torch
import torch.nn.functional as F
from transformers import AutoTokenizer, AutoModelForMaskedLM
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence", required=True, type=str)
    parser.add_argument("--cdr_indices", required=True, type=str)
    parser.add_argument("--num_seqs", required=True, type=int)
    parser.add_argument("--out_fasta", required=True, type=str)
    parser.add_argument("--omit_aas", default="", type=str)
    args = parser.parse_args()

    clean_wt = args.sequence.strip()
    omit_set = set(args.omit_aas.strip().upper())

    # Parse the comma-separated string back into a Python set
    if args.cdr_indices.strip():
        combined_cdr_mask = {int(x) for x in args.cdr_indices.split(",") if x.strip()}
    else:
        combined_cdr_mask = set()

    alphabet_list = list("ACDEFGHIKLMNPQRSTVWY")

    # Load the official Biohub checkpoint (native to modern transformers)
    hf_id = "biohub/ESMC-600M-hf"
    tokenizer = AutoTokenizer.from_pretrained(hf_id)
    model = AutoModelForMaskedLM.from_pretrained(hf_id).eval()
    mask_token_id = tokenizer.mask_token_id

    # Use half-precision on GPU, with a clean CPU fallback if VRAM is exhausted
    device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cuda":
        try:
            amp_dtype = torch.bfloat16 if torch.cuda.is_bf16_supported() else torch.float16
            model = model.to(device=device, dtype=amp_dtype)
        except torch.OutOfMemoryError:
            print(" -> GPU OOM detected. Falling back to CPU execution for ESM-C...")
            torch.cuda.empty_cache()
            device = "cpu"
            model = model.to(device=device, dtype=torch.float32)
    else:
        model = model.to(device=device, dtype=torch.float32)

    inputs = tokenizer(clean_wt, return_tensors="pt")
    inputs = {k: v.to(device) for k, v in inputs.items()}

    tokens = tokenizer.convert_ids_to_tokens(inputs["input_ids"][0])
    special_tokens = set(tokenizer.all_special_tokens or [])
    seq_start_idx = 0
    for idx, tok in enumerate(tokens):
        if tok in special_tokens:
            continue
        if tok.upper() == clean_wt[0].upper():
            seq_start_idx = idx
            break

    plm_muts = []

    for i, wt_aa in enumerate(clean_wt):
        if i in combined_cdr_mask:
            continue

        token_pos = i + seq_start_idx

        masked_input_ids = inputs["input_ids"].clone()
        masked_input_ids[0, token_pos] = mask_token_id

        with torch.no_grad():
            logits = model(input_ids=masked_input_ids).logits[0]

        log_probs = F.log_softmax(logits.float(), dim=-1)

        wt_tok_id = tokenizer.convert_tokens_to_ids(wt_aa)
        wt_prob = log_probs[token_pos, wt_tok_id].item()

        for mut_aa in alphabet_list:
            if mut_aa == wt_aa or mut_aa in omit_set:
                continue
            mut_tok_id = tokenizer.convert_tokens_to_ids(mut_aa)

            if mut_tok_id == tokenizer.unk_token_id:
                continue

            mut_prob = log_probs[token_pos, mut_tok_id].item()
            delta_log_p = mut_prob - wt_prob

            if delta_log_p > 0:
                plm_muts.append({
                    'pos': i, 'wt_aa': wt_aa, 'mut_aa': mut_aa, 'delta': delta_log_p
                })

    plm_muts = sorted(plm_muts, key=lambda x: x['delta'], reverse=True)
    top_muts = plm_muts[:args.num_seqs]

    plm_records = []
    for i, m in enumerate(top_muts):
        mut_seq = list(clean_wt)
        mut_seq[m['pos']] = m['mut_aa']
        new_record = SeqRecord(
            Seq("".join(mut_seq)),
            id=f"ESMC_{m['wt_aa']}{m['pos']+1}{m['mut_aa']}_rank{i+1}",
            description=""
        )
        plm_records.append(new_record)

    if plm_records:
        SeqIO.write(plm_records, args.out_fasta, "fasta")


if __name__ == "__main__":
    main()
