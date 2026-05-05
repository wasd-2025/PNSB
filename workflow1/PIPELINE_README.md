# DSM123 Metabolic Reconstruction Project

This project has been reorganized into a standard Python layout while keeping notebook logic and thresholds unchanged.

## Project structure

```text
D:\Project\ai_project\DSM123\workflow1
├── src\dsm123_pipeline\
│   ├── __init__.py
│   ├── config.py
│   ├── pipeline.py
│   └── cli.py
├── scripts\run_pipeline.py
├── metabolic_network_reconstruction_project.py  (compat entrypoint)
├── genomes\  (put .gb files here)
├── prots\
├── nucl\
├── bbh\
└── Models\
```

## Where to put your new `.gb` file

Recommended location:

- `D:\Project\ai_project\DSM123\workflow1\genomes\<your_target_id>.gb`

Example:

- `D:\Project\ai_project\DSM123\workflow1\genomes\new_strain.gb`

You can also pass any external path with `--target-gb`; the pipeline will copy it into `genomes\<target_id>.gb`.

## Run pipeline

From `D:\Project\ai_project\DSM123\workflow1`:

```bash
python scripts/run_pipeline.py --target-id new_strain --target-gb "D:\path\to\new_strain.gb"
```

or use the root compatible entrypoint:

```bash
python metabolic_network_reconstruction_project.py --target-id new_strain --target-gb "D:\path\to\new_strain.gb"
```

## Overwrite behavior

By default, the run is overwrite mode.

- Existing generated files are removed first.
- Newly generated outputs use the same filenames and overwrite old results.
- This guarantees the files are newly produced for the current `.gb` input.

Disable cleanup only if needed:

```bash
python scripts/run_pipeline.py --target-id new_strain --target-gb "D:\path\to\new_strain.gb" --no-overwrite
```

## Preserved thresholds

- BBH coverage: `COV >= 0.2`
- Orthology binarization: `PID > 65.0`
- BLASTn filter: `PID > 70` and `alnLength > 0.8 * queryEnd`
- Unannotated ORF candidate: `PID >= 80`
- BLAST `evalue`: `0.001`
