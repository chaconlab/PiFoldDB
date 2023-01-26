# PiFoldDB

This is a curated CATH 4.3 dataset for [PiFold](https://github.com/A4Bio/PiFold) (updated version of CATH 4.2 by [Ingraham, et al, NeurIPS 2019](https://github.com/jingraham/neurips19-graph-protein-design)).
This new version included better structures (PDB-REDO), more chains, last CATH release, removed changes with large missing regions (XXX), and included CB (GLY has a fake one). 

## Datasets

Preprocessed data and splits files are:

- [chain_set500.jsonl.gz](chain_set500.jsonl.gz)   Max sequence lenght 500 aa
- [chain_set800.jsonl.gz](chain_set800.jsonl.gz)   Max sequence lenght 800 aa
- [splits_set500.json.gz](splits_set500.json.gz)   Test: 1414 Train: 19164 Validation: 1323
- [splits_set800.json.gz](splits_set800.json.gz)   Test: 1432 Train: 20513 Validation: 1550


