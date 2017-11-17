echo ">>> G-QUADRUPLEXES AROUND 5P <<<"
bin/htsdb-gquad-pred \
	--db Arkon5/db/sqlite.db \
	--table transcr \
	--fasta transcripts.fa \
	--bed cds.regions.on.transcript.bed \
	--pos 5p \
	--span-up 100 \
	--span-down 100 \
	--collapse \
	> gquad_pred.5p.cds.tab

plot-gquad-rel-pos.R \
	--ifile gquad_pred.5p.cds.tab \
	--figfile gquad_pred.5p.cds.pdf
