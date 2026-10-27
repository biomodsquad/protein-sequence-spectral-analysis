# CPP dataset provenance and use

Reviewed 2026-10-27. The maintainer identified the source studies as:

- Manavalan B, Subramaniyam S, Shin TH, Kim MO, Lee G (2018).
  *Machine-Learning-Based Prediction of Cell-Penetrating Peptides and Their
  Uptake Efficiency with Improved Accuracy*. Journal of Proteome Research
  17(8), 2715–2726. [doi:10.1021/acs.jproteome.8b00148](https://doi.org/10.1021/acs.jproteome.8b00148).
- Manavalan B, Patra MC (2022). *MLCPP 2.0: An Updated Cell-penetrating Peptides
  and Their Uptake Efficiency Predictor*. Journal of Molecular Biology
  434(11), 167604. [doi:10.1016/j.jmb.2022.167604](https://doi.org/10.1016/j.jmb.2022.167604).

The active example uses the complete **Layer1** files from the official
[MLCPP 2.0 download page](https://balalab-skku.org/mlcpp2/download/).
The historical repository CSVs are not used. Layer1 labels describe CPP/non-CPP
classification, not Layer2 uptake efficiency or ABP_GN antibacterial activity.

| Official download | Rows | Positive | Negative |
| --- | ---: | ---: | ---: |
| Training | 1,146 | 573 | 573 |
| Independent test | 2,341 | 157 | 2,184 |
| Test after exact training-overlap exclusion | 2,337 | 153 | 2,184 |

## Pinned website snapshot

[mlcpp2-source.json](mlcpp2-source.json) records the official URLs, SHA-256
hashes, class counts, and snapshot date (2026-10-27):

- [Layer1 training FASTA](https://balalab-skku.org/mlcpp2/download/Layer1_training.txt/txt/).
- [Layer1 independent FASTA](https://balalab-skku.org/mlcpp2/download/Layer1_Independent.txt/txt/).

Prepare data directly from the website:

```bash
python research/prepare_cpp_data.py --output-dir research/data/mlcpp2
```

The helper validates **both** downloads against the committed hashes before
writing output. A changed website file fails verification and requires an
explicit snapshot/provenance review. It preserves all official records, labels,
FASTA identifiers, and row order; no negative subsampling, sorting, or length
filtering is applied. Conversion produces training.csv and independent.csv,
plus the original FASTAs and source.json with both raw and prepared-file hashes.

For offline reproduction from the previously verified downloads:

```bash
python research/prepare_cpp_data.py \
  --source-dir .research-deps \
  --output-dir research/data/mlcpp2
```

The source cache must contain mlcpp2_training.fasta and
mlcpp2_independent.fasta with exactly the recorded bytes. The output directory
must be new. Each workflow run also checks the prepared CSV hashes against
source.json and embeds the complete source metadata in report.json.

## Exact overlap handling

Neither official file contains duplicate sequences within itself. Four positive
peptides occur in both training and test files, with matching labels. Training
is preserved; these four test rows are excluded before evaluation. The resulting
test set has 2,337 rows: 153 positive and 2,184 negative. The exclusions record
official FASTA identifiers, original test row indices, and sequence hashes.
All other official test records are evaluated. The downloaded files and
prepared CSVs retain the four overlaps; only the evaluation view excludes them.

Unknown/empty sequences, repeated rows within a file, duplicate official
identifiers, or conflicting cross-file labels raise errors rather than being
silently imputed or deduplicated. Prediction and split manifests retain official
identifiers. There is no family-disjointness guarantee, and exact-overlap removal
does not make this a new independent biological validation. Scores are not
identical to a paper evaluation that retains a different test protocol.

## Data distribution

Raw downloads and prepared CSVs are local cache files in ignored research/data/.
They are excluded from both wheel and source distributions; source URLs, hashes,
preparation code, and instructions are included. Reports contain identifiers
rather than additional raw peptide sequences.

The official download page does not state an explicit redistribution license in
the page inspected on the review date. The code's BSD license is not a statement
of third-party dataset rights. This example obtains the data from its source
for analysis and makes no separate dataset redistribution claim.
