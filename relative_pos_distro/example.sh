echo ">>>> RELATIVE POS DISTRO 5p-5p <<<<"
bin/relative-pos-distro-on-feats \
	--db1 Akron5/db/sqlite.db \
	--table1 transcr \
	--pos1 5p \
	--collapse1 \
	--db2 RiboSeq/db/sqlite.db \
	--table2 transcr \
	--pos2 5p \
	--collapse2 \
	--bed cds.regions.on.transcript.bed \
	--span 60 \
	--offset 100 \
	> rel-pos.5p_5p.cds.tab

src/R/rel-pos-distro.R \
	--ifile rel-pos.5p_5p.cds.tab \
	--figfile rel-pos.5p_5p.cds.pdf

echo ">>>> RELATIVE POS DISTRO 5p-5p PER TRANSCRIPT <<<<"
bin/relative-pos-distro-on-feats \
	--db1 Akron5/db/sqlite.db \
	--table1 transcr \
	--pos1 5p \
	--collapse1 \
	--db2 RiboSeq/db/sqlite.db \
	--table2 transcr \
	--pos2 5p \
	--collapse2 \
	--bed cds.regions.on.transcript.bed \
	--span 60 \
	--offset 100 \
	--by-feat \
	> rel-pos.5p_5p.cds.per_transcr.tab
bin/relative-pos-distro-on-feats \
	--db1 Akron5/db/sqlite.db \
	--table1 transcr \
	--pos1 5p \
	--collapse1 \
	--random1 \
	--db2 RiboSeq/db/sqlite.db \
	--table2 transcr \
	--pos2 5p \
	--collapse2 \
	--bed cds.regions.on.transcript.bed \
	--span 60 \
	--offset 100 \
	--by-feat \
	> rel-pos.rand5p_5p.cds.per_transcr.tab

echo ">>>> MAKE HEATMAP <<<<"
rel-pos-distro-heatmap.R \
	--ifile rel-pos.5p_5p.cds.per_transcr.tab \
	--rfile rel-pos.rand5p_5p.cds.per_transcr.tab \
	--quantiles 10 \
	--posMin 40 \
	--posMax 40 \
	--color steelblue \
	--figfile rel-pos.5p_5p.cds.heatmap10.pdf
