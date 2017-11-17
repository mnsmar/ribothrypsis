package main

import (
	"bytes"
	"fmt"
	"io"
	"log"
	"math"
	"math/rand"
	"os"
	"regexp"

	_ "github.com/mattn/go-sqlite3"

	"github.com/Masterminds/squirrel"
	"github.com/alexflint/go-arg"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/io/featio"
	"github.com/biogo/biogo/io/featio/bed"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/jmoiron/sqlx"
)

// Opts is the struct with the options that the program accepts.
type Opts struct {
	DB       string `arg:"required,help:SQLite3 database"`
	Table    string `arg:"required,help:database table name"`
	Where    string `arg:"help:SQL filter injected in WHERE clause"`
	FASTA    string `arg:"required,help:FASTA file with reference sequences"`
	BED      string `arg:"required,help:BED file with features"`
	Pos      string `arg:"required,help:reference point for reads; one of 5p or 3p"`
	SpanUp   int    `arg:"--span-up,required,help:max upstream distance from pos"`
	SpanDown int    `arg:"--span-down,required,help:max downstream distance from pos"`
	Collapse bool   `arg:"help:collapse reads that have the same pos."`
	LenLim   int    `arg:"--len-lim,help:minimum feature length"`
	QArea    int    `arg:"--q-area,help:area (nts) to extend query around feats (Default: 1000)"`
	Random   bool   `arg:"help:randomize pos"`
	Shuffle  bool   `arg:"help:shuffle sequence around pos maintaining codon composition"`
	Orf      int    `arg:"help:keep G-quadruplexes on this ORF only; disabled when negative"`
}

// Version returns the program version.
func (Opts) Version() string {
	return `htsdb-gquad-pred 0.1`
}

// Description returns an extended description of the program.
func (Opts) Description() string {
	return `Prints the distribution of G-Quadruplexes around the 5'/3' end of reads that are located in features defined in a BED file. G-Quadruplex positions are predicted with the pattern G{2,4}.{1,7}G{2,4}.{1,7}G{2,4}.{1,7}G{2,4}`
}

func main() {
	// Define command line options.
	var opts Opts
	opts.QArea = 1000
	opts.Orf = -1
	p := arg.MustParse(&opts)
	if opts.Pos != "5p" && opts.Pos != "3p" {
		p.Fail("--pos must be either 5p or 3p")
	}
	if opts.SpanUp < 0 || opts.SpanDown < 0 {
		p.Fail("--span-up and --span-down must be positive integers")
	}

	// Define the G-Quadruplex prediction pattern.
	validG4 := regexp.MustCompile(`G{2,4}.{1,7}G{2,4}.{1,7}G{2,4}.{1,7}G{2,4}`)
	margin := 50 // min sequence length upstream and downstream

	// Get position extracting function.
	extractPos := Head
	if opts.Pos == "3p" {
		extractPos = Tail
	}

	// Read FASTA and store sequences in a map.
	fastaF, err := os.Open(opts.FASTA)
	if err != nil {
		log.Fatal(err)
	}
	defer Close(fastaF)
	fastaR := fasta.NewReader(fastaF, linear.NewSeq("", nil, alphabet.DNA))
	fastaS := seqio.NewScanner(fastaR)
	seqs := make(map[string][]byte)
	for fastaS.Next() {
		s := fastaS.Seq().(*linear.Seq)
		l := alphabet.LettersToBytes(s.Seq)
		seqs[s.Name()] = bytes.ToUpper(l)
	}
	if err = fastaS.Error(); err != nil {
		log.Fatal(err)
	}

	// Open BED6 scanner
	bedF, err := os.Open(opts.BED)
	if err != nil {
		log.Fatal(err)
	}
	defer Close(bedF)
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

	// Prepare squirrel select.
	sel := squirrel.Select("strand", "rname", "start", "stop").From(opts.Table).Where(
		"strand = ? AND rname = ? AND start BETWEEN ? AND ? AND stop BETWEEN ? AND ?")
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

	// Loop on the features.
	counts := make(map[int]int)
	recordCnt := 0
	for bedS.Next() {
		b := bedS.Feat()

		if opts.LenLim != 0 && b.Len() < opts.LenLim {
			continue
		}

		o, ok := b.(feat.Orienter)
		if !ok {
			log.Fatal("error: no orientation for bed feature")
		}
		if o.Orientation() == feat.Reverse {
			log.Fatal("error: bed features on reverse orientation are not supported")
		}

		qOri := o.Orientation()
		qStart := b.Start() - opts.QArea
		qStop := b.End() - 1 + opts.QArea
		qRname := b.Location().Name()
		bed5p := b.Start()
		bed3p := b.End() - 1

		seq, ok := seqs[qRname]
		if !ok {
			continue
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

			// pos with span should be within bed5p and bed3p
			if pos < bed5p+opts.SpanUp+margin || pos > bed3p-opts.SpanDown-margin {
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
			Randomize(wig, bed5p+opts.SpanUp+margin, bed3p-opts.SpanDown-margin+1)
		}

		for pos, val := range wig {
			from := pos - opts.SpanUp - margin
			to := pos + opts.SpanDown + margin + 1
			subseq := seq[from:to]
			if opts.Shuffle {
				var orfOffset int
				switch math.Mod(float64(from-bed5p), 3) {
				case 0:
					orfOffset = 0
				case 1:
					orfOffset = 2
				case 2:
					orfOffset = 1
				}
				for i := orfOffset; i < len(subseq)-2; i += 3 {
					var j int
					for {
						j = rand.Intn(i + 1)
						if math.Mod(float64(i-j), 3) == 0 {
							break
						}
					}
					subseq[i], subseq[j] = subseq[j], subseq[i]
					subseq[i+1], subseq[j+1] = subseq[j+1], subseq[i+1]
					subseq[i+2], subseq[j+2] = subseq[j+2], subseq[i+2]
				}
			}
			matches := validG4.FindAllIndex(subseq, -1)
			for _, m := range matches {
				matchStart, matchEnd := m[0], m[1]
				matchPos := from + matchStart
				orf := (matchPos - bed5p) % 3
				if opts.Orf > -1 && opts.Orf != orf {
					continue
				}
				for i := matchStart; i < matchEnd; i++ {
					counts[i-opts.SpanUp-margin] += val
				}
			}
			recordCnt += val
		}
	}
	if err = bedS.Error(); err != nil {
		log.Fatal(err)
	}

	// Print the counts.
	fmt.Print("pos\tval\trecords\n")
	for pos := -opts.SpanUp; pos <= opts.SpanDown; pos++ {
		fmt.Printf("%d\t%d\t%d\n", pos, counts[pos], recordCnt)
	}
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

// Close calls the Close() for a Closer and checks the returned error. If
// an error is returned it is logged and os.Exit(1) is called.
func Close(c io.Closer) {
	if err := c.Close(); err != nil {
		log.Fatal(err)
	}
}

// Randomize randomizes wig using Sattolo's algorithm. Affected keys will be
// in the region [from, to).
func Randomize(wig map[int]int, from, to int) {
	for i := from; i < to; i++ {
		j := from + rand.Intn(i-from+1)
		wig[i], wig[j] = wig[j], wig[i]
	}
}
