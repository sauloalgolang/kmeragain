package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"sort"
	"strings"
)

const INVB byte = '\x7F'
const INVR rune = rune(INVB)

type LISTADDSTATUS byte

const LISTADDOK LISTADDSTATUS = 0
const LISTADDFULL LISTADDSTATUS = 1

type Unit struct {
	seq   uint64
	count uint16
}

type unitSlice []Unit
type UnitList struct {
	units    unitSlice
	count    uint64
	total    uint64
	min      uint64
	max      uint64
	blockMin uint64
	blockMax uint64
}

type unitGroupPointer []*UnitList
type unitGroup []UnitList
type UnitListDb struct {
	lists     unitGroupPointer
	repo      unitGroup
	count     uint64
	min       uint64
	max       uint64
	blockMin  uint64
	blockMax  uint64
	blockSize uint64 //TODO: make dynamic
	blockNum  uint64
}

func (units *unitSlice) Append(unit Unit) {
	(*units) = append((*units), unit)
}

func (units *unitSlice) Sort() {
	sort.Slice((*units)[:], func(i, j int) bool {
		return (*units)[i].seq < (*units)[j].seq
	})
}

func (list *UnitList) Add(seq uint64) LISTADDSTATUS {
	if seq >= list.min {
		for pos, unit := range list.units {
			// TODO: BINARY SEARCH
			if seq == unit.seq {
				list.units[pos].count += 1
				list.count += 1
				return LISTADDOK
			} else if seq > unit.seq {
				break
			}
		}
	}

	if seq < list.min {
		list.min = seq
	}

	if seq > list.max {
		list.max = seq
	}

	list.count += 1
	list.total += 1
	list.units.Append(Unit{seq, 1})
	list.units.Sort()

	// TODO: SPLIT IF TOO LARGE
	if list.blockMin > list.blockMax {
		return LISTADDFULL
	} else {
		return LISTADDOK
	}
}

func (list *unitGroupPointer) Sort() {
	sort.Slice((*list)[:], func(i, j int) bool {
		return (*list)[i].blockMin < (*list)[j].blockMin
	})
}

func (list *unitGroupPointer) Append(repo *unitGroup) {
	*list = append(*list, &(*repo)[len(*repo)-1])
}

func (db *UnitListDb) Create(seq uint64) {
	if seq > db.max {
		db.max = seq
	}

	if seq < db.min {
		db.min = seq
	}

	// blockSize uint64 //TODO: make dynamic
	// blockNum  uint64

	// TODO: CREATE
	//  Calculate blocks
	//    update list
	//    update db
	// TODO: SPLIT IF TOO LARGE

	newul := UnitList{}
	db.repo = append(db.repo, newul)
	db.lists.Append(&db.repo)
	db.lists.Sort()
}

func (db *UnitListDb) Add(seq uint64) {
	if seq < db.blockMin || seq > db.blockMax {
		db.Create(seq)

	} else {
		for pos, list := range db.lists {
			// TODO: BINARY SEARCH
			if seq >= list.blockMin {
				if seq <= list.blockMax {
					db.lists[pos].Add(seq)
					return
				} else {
					break
				}
			}
		}
		db.Create(seq)
	}
}

func main() {
	filename := "S_lycopersicum_chromosomes.2.50.fa"
	// filename := "test.fasta"
	fragsize := 21

	file, err := os.Open(filename)

	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	fmt.Println("reading")

	scanner := bufio.NewScanner(file)

	// leftoverLine := ""
	leftoverDigi := ""
	row := 0
	valids := 0
	invalids := 0
	seqnum := -1

	var countvalids []int
	var countinvalids []int
	var names []string
	var line string
	var digi string
	var linelen int
	var linepos int
	var res int
	var valid bool
	for scanner.Scan() {
		//fmt.Println("leftover 1", leftover)

		line = scanner.Text()

		//fmt.Println(line) // Println will add back the final '\n'
		//fmt.Println("leftover 2", leftover)

		row += 1

		if len(line) == 0 {
			leftoverDigi = ""
		} else if line[0] == '>' {
			//leftoverLine = ""
			leftoverDigi = ""
			countvalids = append(countvalids, 0)
			countinvalids = append(countinvalids, 0)
			names = append(names, line[1:])
			seqnum += 1
			fmt.Println(line)
		} else {
			//fmt.Println("leftover before", leftover)
			//fmt.Println("line     before", line, len(line))
			//line =          strings.ToUpper(line)
			digi = dna2int(strings.ToUpper(line))
			//line = leftoverLine + line
			digi = leftoverDigi + digi
			//fmt.Println("line     after ", line, len(line))
			linelen = len(digi)

			for linepos = 0; linepos < linelen-fragsize+1; linepos++ {
				fragDigi := digi[linepos : linepos+fragsize]

				/*
					fragLine := line[linepos:linepos+fragsize]

					fmt.Println(row, linepos, len(digi))

					for i := 0; i < len(line); i++ {
							fmt.Printf("%s ", string(line[i]))
						}
					fmt.Print(" ")
					for i := 0; i < len(fragLine); i++ {
							fmt.Printf("%s ", string(fragLine[i]))
						}
					fmt.Println("")

					for i := 0; i < len(digi); i++ {
							fmt.Printf("%x ", digi[i])
						}
					fmt.Print(" ")
					for i := 0; i < len(fragDigi); i++ {
							fmt.Printf("%x ", fragDigi[i])
						}
					fmt.Println("\n")

					fragLine = fragLine + ""
				*/

				valid, res = frag2num(fragDigi)

				// fmt.Printf("% 07x %t %d %014b\n", fragDigi, valid, res, res)

				if valid {
					res += 0
					valids += 1
					countvalids[seqnum] += 1
				} else {
					invalids += 1
					countinvalids[seqnum] += 1
				}
			}

			//leftoverLine = string(line[linelen-fragsize+1:])
			leftoverDigi = string(digi[linelen-fragsize+1:])
			//fmt.Println("leftover after ", ""+leftover+"")

			// fmt.Println("")
		}
	}

	if err := scanner.Err(); err != nil {
		fmt.Fprintln(os.Stderr, "error reading file:", err)
	}

	totalvalids := 0
	totalinvalids := 0
	numseqs := 0
	for i := 0; i < len(names); i++ {
		fmt.Printf("%-20s %12d %12d\n", names[i], countvalids[i], countinvalids[i])
		numseqs += 1
		totalvalids += countvalids[i]
		totalinvalids += countinvalids[i]
	}
	fmt.Println("")
	fmt.Printf("%-20s %12d\n", "seqs", numseqs)
	fmt.Printf("%-20s %12d %12d %12d\n", "total", totalvalids, totalinvalids, totalvalids+totalinvalids)
	fmt.Printf("%-20s %12d %12d %12d\n", "count", valids, invalids, valids+invalids)
}

// https://gobyexample.com/collection-functions
func Map(vs string, f func(rune) rune) string {
	vsm := make([]rune, len(vs))
	for i, v := range vs {
		vsm[i] = f(v)
	}
	return string(vsm)
}

func dna2int(line string) string {
	return Map(line, nuc2int)
}

func nuc2int(nuc rune) rune {
	switch nuc {
	case 'A':
		return '\x00'
	case 'C':
		return '\x01'
	case 'G':
		return '\x02'
	case 'T':
		return '\x03'
	default:
		return INVR
	}
}

func frag2num(fragDigi string) (bool, int) {
	/*
		TODO:
		- return false if has value larger then 3
		- Get 4 by 4 from left to right
		- Get fwd and rev value for byte
		- assign value to correct position in fwd and rev arrays
		- return smallest value
	*/

	/*
		for i := 0; i < len(fragDigi); i++ {
			fmt.Printf("%x ", fragDigi[i])
		}
	*/

	val := 0

	// if strings.IndexByte(fragDigi, INVB) != -1 {
	// 	return false, 0
	// } else {
	ld := len(fragDigi)
	for i := 0; i < ld; i++ {
		digi := fragDigi[i]

		if digi == INVB {
			return false, 0
		}

		nval := int(digi) << (uint(ld-i-1) * 2)

		val += nval

		// fmt.Printf(" %d %02x %012d %014b\n", i, digi, nval, nval)
	}

	// fmt.Printf("      %012d %014b\n", val, val)

	return true, val
	// }
}
