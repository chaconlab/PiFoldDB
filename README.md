# PiFoldDB

This is a curated CATH 4.3 dataset for [PiFold](https://github.com/A4Bio/PiFold) (an updated version of CATH 4.2 by [Ingraham, et al, NeurIPS 2019](https://github.com/jingraham/neurips19-graph-protein-design)).
This new version included better structures (PDB-REDO), more chains, the last CATH release, included gaps (noted by "-") and missing regions (noted as "X" with NaN coordinates), removed tags, and cases with large missing regions.  

## Datasets

Preprocessed data and splits can be found here: [cathPi.tgz](https://saco.csic.es/index.php/s/JR5BaqeiBGg7G3D):

- chain_set.jsonl   Max sequence length 500 aa
- chain_set_splits.json   Test: 1422 Train: 18960 Validation: 1436




