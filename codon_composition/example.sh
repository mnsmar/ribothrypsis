echo ">>>> CODON DISTRIBUTION AROUND 5p <<<<"
codon_composition \
	--db Akron5/db/sqlite.db \
	--table transcr \
	--fasta transcripts.fa \
	--collapse \
	--pos 5p \
	--span-up 20 \
	--span-down 20 \
	--bed cds.regions.on.transcript.bed \
| add_aminoacid.pl \
	--ifile - \
	--codon-col codon \
	--out-col aa \
| table-paste-col.pl --table - --col-name cond --col-val real \
> codon_composition.5p.tab

echo ">>>> CODON DISTRIBUTION AROUND 5p (RANDOM - KEEP CONTENT) <<<<"
codon_composition \
	--db Akron5/db/sqlite.db \
	--table transcr \
	--fasta transcripts.fa \
	--collapse \
	--pos 5p \
	--span-up 20 \
	--span-down 20 \
	--win-up 1 \
	--win-down 0 \
	--random \
	--bed cds.regions.on.transcript.bed \
| add_aminoacid.pl \
	--ifile - \
	--codon-col codon \
	--out-col aa \
| table-paste-col.pl --table - --col-name cond --col-val random \
> codon_composition.5p.random.tab

echo ">>>> MAKE PLOT TO COMPARE REAL TO RANDOM <<<<"
table-cat.pl \
	codon_composition.5p.tab \
	codon_composition.5p.random.tab \
| codon-enrichment-vs-random.R \
	--ifile stdin \
	--figfile codon_composition.5p_vs_rand5p.pdf
