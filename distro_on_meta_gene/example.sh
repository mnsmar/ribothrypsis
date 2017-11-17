echo ">>>> DISTRIBUTION OF 5p ON META-TRANSCRIPT <<<<"
bin/distro-on-meta-feats \
	--db Akron5/db/sqlite.db \
	--table transcr \
	--bed cds.regions.on.transcript.bed \
	--pos 5p \
	--bins 20 \
	--min-len 200 \
	--collapse \
	--by-feat \
> distro-on-meta-elements.5p.tab

plot-distro-on-meta-elements.R \
	-i distro-on-meta-elements.5p.tab \
	-f distro-on-meta-elements.5p.pdf
