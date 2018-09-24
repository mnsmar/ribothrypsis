* `akron5.test.sqlite.db.gz`: A zipped SQLite3 database with a subset of the AKRON 5 reads. The alignment has been performed on the human transcriptome. To build the database file from a SAM file the command `clipseqtools-preprocess sam_to_sqlite` from clipseqtools can be used. To create the annotation columns `utr5`, `cds`, `utr3` (optional columns) the command `clipseqtools-preprocess annotate_with_file` can be used (the annotation file should be a BED with the coordinates of the corresponding elements on the **transcriptome** (e.g. `cds.regions.on.transcript.bed.gz`).

* `cds.regions.on.transcript.bed.gz`: Annotation file (BED format) with the coordinates of the coding sequence on each transcript.

* `transcripts.fa`: FASTA file with the transcript sequences.
