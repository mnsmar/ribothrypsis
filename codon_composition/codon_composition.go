package main

import (
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"sort"
	"strings"

	_ "github.com/mattn/go-sqlite3"

	"github.com/Masterminds/squirrel"
	"github.com/alexflint/go-arg"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/bed"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/seq/linear"
	"github.com/jmoiron/sqlx"
)

// Opts is the struct with the options that the program accepts.
type Opts struct {
	DB       string `arg:"required,help:SQLite3 database"`
	Table    string `arg:"required,help:database table name"`
	Where    string `arg:"help:SQL filter injected in WHERE clause"`
	FASTA    string `arg:"required,help:FASTA file with reference sequences"`
	BED      string `arg:"required,help:BED file with CDS coordinates (defines features)"`
	Pos      string `arg:"required,help:reference point for reads; one of 5p or 3p"`
	SpanUp   int    `arg:"--span-up,required,help:max upstream distance from pos (in codons)"`
	SpanDown int    `arg:"--span-down,required,help:max downstream distance from pos (in codons)"`
	Collapse bool   `arg:"help:collapse reads that have the same pos."`
	LenLim   int    `arg:"--len-lim,help:minimum feature length"`
	QArea    int    `arg:"--q-area,help:area (nts) to extend query around feats; use a big number Default: 1000"`
	Random   bool   `arg:"help:change to random pos within same features that maintain nucleotide content"`
	WinUp    int    `arg:"--win-up,help:when --random; max upstream distance from pos to maintain nt content (nts)"`
	WinDown  int    `arg:"--win-down,help:when --random; max downstream distance from pos to maintain nt content (nts)"`
	Anti     bool   `arg:"help:select reads on the reverse feature orientation"`
}

// Version returns the program version.
func (Opts) Version() string {
	return "codon-composition 0.3"
}

// Description returns an extended description of the program.
func (Opts) Description() string {
	return "Prints the codon composition of the reference sequence around reads."
}

func main() {
	var opts Opts
	opts.QArea = 1000
	p := arg.MustParse(&opts)
	if opts.Pos != "5p" && opts.Pos != "3p" {
		p.Fail("--pos must be either 5p or 3p")
	}
	if opts.SpanUp < 0 && opts.SpanDown < 0 {
		p.Fail("--span-up and --span-down must be positive integers")
	}
	if opts.WinUp < 0 && opts.WinDown < 0 {
		p.Fail("--win-up and --win-down must be positive integers")
	}

	// get position extracting function
	extractPos := Head
	if opts.Pos == "3p" {
		extractPos = Tail
	}

	// Open FASTA scanner
	fastaS, err := FastaScanner(opts.FASTA)
	if err != nil {
		panic(err)
	}

	// Create a map to store the FASTA seqs.
	seqs := make(map[string]seq.Sequence)
	for fastaS.Next() {
		s := fastaS.Seq()
		seqs[s.Name()] = s
	}

	// open BED6 scanner
	bedF, err := os.Open(opts.BED)
	if err != nil {
		log.Fatal(err)
	}
	bedR, err := bed.NewReader(bedF, 6)
	if err != nil {
		log.Fatal(err)
	}
	bedS := featio.NewScanner(bedR)

	// Connect to the database.
	db, err := sqlx.Connect("sqlite3", opts.DB)
	if err != nil {
		log.Fatal(err)
	}

	// Prepare squirrel select builder.
	sel := squirrel.Select("strand", "rname", "start", "stop").From(opts.Table).Where(
		"strand = ? AND rname = ? AND start BETWEEN ? AND ? AND stop BETWEEN ? AND ?")

	// Apply extra where clause if provided.
	if opts.Where != "" {
		sel = sel.Where(opts.Where)
	}

	// Prepare SQL query.
	query, _, err := sel.ToSql()
	if err != nil {
		log.Fatal(err)
	}

	// Prepare statement.
	stmt, err := db.Preparex(query)
	if err != nil {
		log.Fatal(err)
	}

	// Initialize the counts map.
	counts := make(map[int]map[string]int)
	for i := -opts.SpanUp; i <= opts.SpanDown; i++ {
		counts[i] = make(map[string]int)
	}

	// Loop on the features.
	for bedS.Next() {
		b := bedS.Feat()

		if opts.LenLim != 0 && b.Len() < opts.LenLim {
			continue
		}

		o, ok := b.(feat.Orienter)
		if !ok {
			log.Fatal("error: undefined orientation for feature")
		}
		if o.Orientation() == feat.Reverse {
			log.Fatal("error: features on reverse orientation are unsupported")
		}

		qOri := o.Orientation()
		qStart := b.Start() - opts.QArea
		qStop := b.End() - 1 + opts.QArea
		qRname := b.Location().Name()
		cds5p := b.Start()
		cds3p := b.End() - 1

		seq, ok := seqs[qRname]
		if !ok {
			continue
		}

		if opts.Anti {
			qOri = -1 * qOri
		}
		rows, err := stmt.Queryx(qOri, qRname, qStart, qStop, qStart, qStop)
		if err != nil {
			log.Fatal(err)
		}

		toCollapse := make(map[int]bool)
		wig := make(map[int]int)
		r := &Feature{}
		for rows.Next() {
			if err = rows.StructScan(r); err != nil {
				log.Fatal(err)
			}

			ori := r.Orientation()
			pos := extractPos(r, ori)

			// TODO handle reverse orientation
			// if ori == feat.Reverse {
			// 	vals = reverse(vals)
			// 	lenVals := len(vals)
			// 	pos = lenVals - pos - 1
			// 	cds5p, cds3p = lenVals-cds3p-1, lenVals-cds5p-1
			// }

			// pos with span should be within cds5p and cds3p
			if pos-(3*opts.SpanUp)-2 < cds5p || pos+(3*opts.SpanDown) > cds3p {
				continue
			}

			if opts.Collapse {
				if toCollapse[pos] {
					continue
				}
				toCollapse[pos] = true
			}
			wig[pos]++
		}

		if opts.Random {
			wig = RandomWig(wig, seq, opts.SpanUp, opts.SpanDown, cds5p, cds3p, opts.WinUp, opts.WinDown)
		}

		for pos, val := range wig {
			// pos is the position of the first nucleotide of the codon we are in
			pos = pos - int(math.Mod(float64(pos-cds5p), 3))

			for i := -opts.SpanUp; i <= opts.SpanDown; i++ {
				codonstart := pos + (i * 3)

				if codonstart < cds5p || codonstart > cds3p {
					log.Fatalf("%s: %d", "span out of bounds for pos", pos)
				}

				codon := strings.ToUpper(string(seq.At(codonstart).L))
				codon += strings.ToUpper(string(seq.At(codonstart + 1).L))
				codon += strings.ToUpper(string(seq.At(codonstart + 2).L))

				counts[i][codon] += val
			}
		}
	}
	if err = bedS.Error(); err != nil {
		log.Fatal(err)
	}

	// Print the counts.
	letters := make(map[string]bool)
	for pos := -opts.SpanUp; pos <= opts.SpanDown; pos++ {
		for l := range counts[pos] {
			letters[l] = true
		}
	}
	var lettersSlice []string
	for l := range letters {
		lettersSlice = append(lettersSlice, l)
	}
	sort.Strings(lettersSlice)

	fmt.Print("pos\tcodon\tcount\ttotal_count\n")
	for pos := -opts.SpanUp; pos <= opts.SpanDown; pos++ {
		v := counts[pos]
		totalCnt := 0
		for _, cnt := range v {
			totalCnt += cnt
		}
		for _, l := range lettersSlice {
			fmt.Printf("%d\t%s\t%d\t%d\n", pos, l, v[l], totalCnt)
		}
	}
}

// offPos
type orfPos struct {
	ORF int
	Pos int
}

// Feature is part of an htsdb record that wraps Range and the name of the
// reference.
type Feature struct {
	Rname    string           `db:"rname"`
	Orient   feat.Orientation `db:"strand"`
	StartPos int              `db:"start"`
	StopPos  int              `db:"stop"`
}

// Name returns an empty string.
func (e *Feature) Name() string { return "" }

// Start returns the start position of Range.
func (e *Feature) Start() int { return e.StartPos }

// End returns the end position of Feature.
func (e *Feature) End() int { return e.StopPos + 1 }

// Len returns the length of Feature.
func (e *Feature) Len() int { return e.End() - e.Start() }

// Description returns an empty string.
func (e *Feature) Description() string { return "" }

// Location returns the location of Feature.
func (e *Feature) Location() feat.Feature { return nil }

// Orientation returns the orientation of OrientedFeature.
func (e *Feature) Orientation() feat.Orientation { return e.Orient }

// Head returns the head coordinate of r depending on orientation.
func Head(r feat.Range, o feat.Orientation) int {
	if o == feat.Forward {
		return r.Start()
	} else if o == feat.Reverse {
		return r.End() - 1
	}
	panic("orientation must be forward or reverse")
}

// Tail returns the tail coordinate of r depending on orientation.
func Tail(r feat.Range, o feat.Orientation) int {
	if o == feat.Forward {
		return r.End() - 1
	} else if o == feat.Reverse {
		return r.Start()
	}
	panic("orientation must be forward or reverse")
}

// FastaScanner returns a seqio.Scanner that reads from f.
func FastaScanner(f string) (*seqio.Scanner, error) {
	ioR, err := os.Open(f)
	if err != nil {
		return nil, err
	}
	fR := fasta.NewReader(ioR, linear.NewSeq("", nil, alphabet.DNA))
	return seqio.NewScanner(fR), nil
}

// RandomWig returns a randomized wig keeping the same nucleotide content.
func RandomWig(
	wig map[int]int, sequence seq.Sequence, SpanUp, SpanDown, cds5p, cds3p, WinUp, WinDown int) map[int]int {

	var key string
	rwig := make(map[int]int)
	lseq := strings.ToUpper(sequence.(*linear.Seq).String())

	// Initialize the counts map.
	content := make(map[orfPos][]string)

	seqCnt := make(map[int]int)
	for p, v := range wig {
		orf := int(math.Mod(float64(p-cds5p), 3))
		seqCnt[orf] += v
		for i := -WinUp; i <= WinDown; i++ {
			pos := p + i
			nt := string(lseq[pos])
			s := content[orfPos{ORF: orf, Pos: i}]
			for j := 0; j < v; j++ {
				s = append(s, nt)
			}
			content[orfPos{ORF: orf, Pos: i}] = s
		}
	}

	index := make(map[string]map[int][]int)
	for p := cds5p + 3*SpanUp + 2; p <= cds3p-3*SpanDown; p++ {
		orf := int(math.Mod(float64(p-cds5p), 3))
		key = ""
		for i := -WinUp; i <= WinDown; i++ {
			pos := p + i
			key += string(lseq[pos])
		}
		if _, ok := index[key]; !ok {
			index[key] = make(map[int][]int)
		}
		index[key][orf] = append(index[key][orf], p)
	}

	var tries int
	for orf := 0; orf < 3; orf++ {
		for seqCnt[orf] > 0 && tries < 1000 {
			key = ""
			for i := -WinUp; i <= WinDown; i++ {
				s := content[orfPos{ORF: orf, Pos: i}]
				key += string(s[rand.Intn(len(s))])
			}

			positions, ok := index[key][orf]
			if !ok {
				tries++
				continue
			}
			tries = 0
			seqCnt[orf]--

			randomPos := positions[rand.Intn(len(positions))]
			rwig[randomPos]++
		}
	}

	return rwig
}
