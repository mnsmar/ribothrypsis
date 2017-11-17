echo ">>>> RELATIVE POS 5p-STOP CODON <<<<"
bin/dist-from-feats \
	--db Akron5/db/sqlite.db \
	--table transcr \
	--bed cds.regions.on.transcript.bed \
	--bed-pos 3p \
	--pos 5p \
	--span-up 99 \
	--span-down 99 \
	--collapse \
	> rel-pos-5p-stop.tab

plot-dist.R \
	--ifile rel-pos-5p-stop.tab \
	--figfile rel-pos-5p-stop.pdf \
	--posMin '-97' --posMax '99' --riboPos '-17' \
	--xPosMin '-80'	--xPosMax '20'
