# alignment_to_bam, Explained from Scratch

## Part 1: The Biology

### DNA is a very long text
Every living thing carries instructions called **DNA**. You can picture DNA as a very long book written with only four letters:

**A   C   G   T**

A human's DNA is about **3 billion letters** long, split into 23 or so "chapters" called **chromosomes** (chr1, chr2, …).

### DNA has two sides, called strands
DNA is like a zipper with two sides. The letters on one side determine the other side, because they always pair up:

- **A** pairs with **T**
- **C** pairs with **G**

So if one side reads `ACCG`, the other side is `TGGC`, and it runs in the opposite direction. Read in its own direction, the other side says **`CGGT`**.

That flipped version is called the **reverse complement**: you reverse the order, then swap each letter for its partner. We'll need this later.

### Machines can only read small pieces
We can't read all 3 billion letters in one go. A **sequencing machine** shreds the DNA into millions of tiny pieces and reads each one, about 100–150 letters at a time. Each piece is called a **read**.

It's like shredding a book into millions of strips of paper and reading each strip.

Along with every letter, the machine records how sure it was about that letter. These are **quality scores**.

### Paired reads
Often the machine reads **both ends** of the same small piece of DNA. That gives two reads, called **mates**, which belong together. They're usually named something like `read42/1` and `read42/2`.

### Alignment: putting the strips back in place
We already have a "master copy" of the human book, called the **reference genome**. For each strip, we want to answer:

> "Where in the master book does this strip come from?"

Answering that is called **aligning** or **mapping** the read. A complete answer says:

- **which chapter** (chromosome), e.g. `chr1`;
- **which position**, e.g. letter number 1,234,567;
- **which side of the zipper** (strand): forward or reverse;
- **how well it fits**. Strips rarely match perfectly, because people differ slightly and machines make mistakes. So we also record where letters match, where there are extra letters, and where letters are missing.

That last bit gets written in a compact code called **CIGAR**:

`50M 2I 48M` = "50 letters match, 2 extra letters, then 48 more match"

| CIGAR letter | Meaning |
|---|---|
| M | The letters line up (match or mismatch) |
| I | Insertion: the read has extra letters the book doesn't have |
| D | Deletion: the read is missing letters the book has |
| S | Soft clip: letters at the edge of the read that didn't fit anywhere |

## Part 2: What vg Is and Why This Function Exists

Most tools compare strips to **one** master book. **vg** is smarter: it uses a **graph genome**. Picture a book where some sentences have several versions printed side by side, covering the common ways different people's DNA differs. That makes it better at placing strips.

The downside is that the rest of the science world uses the plain one-book system. Everyone shares results in a standard file format called **SAM** (readable text) or **BAM** (the same thing, compressed into binary, which takes less space and is faster).

So vg's job, at the end, is to **translate its answers into the format everyone else understands**.

That's what `alignment_to_bam` does. It takes vg's answer for **one read** and fills in **one standard BAM record** for it.

> **Analogy:** You've written personal notes about a package ("blue box, going to Ana, 3rd floor"). The post office only accepts its own official form, with fixed boxes in a fixed order. `alignment_to_bam` copies your notes onto the official form, in exactly the right boxes, in exactly the right format.

The function **doesn't figure out where the read goes**. Another part of vg already did that. This function only **fills in the form**.

## Part 3: The Inputs, in Plain Language

```
bam1_t* alignment_to_bam(bam_hdr_t* bam_header,
                         const Alignment& alignment,
                         const string& refseq,
                         const int32_t refpos,
                         const bool refrev,
                         const vector<pair<int, char>>& cigar);
```

| Input | Plain meaning | Example |
|---|---|---|
| bam_header | The file's table of contents: a list of chapter names | chr1, chr2, … |
| alignment | vg's notes about the read: its name, letters, quality scores, score | read42/1, ACGT… |
| refseq | Which chapter it was placed in | "chr1", or "" if it couldn't be placed |
| refpos | Which letter position in that chapter | 1234566, or -1 if not placed |
| refrev | Did it match the other side of the zipper? | true / false |
| cigar | The match description | 50M 2I 48M |

**Output:** a finished BAM record, the completed post-office form.

### A little C++ to read that line
- **string**: text, like "chr1".
- **int32_t**: a whole number (the "32" is its size in memory). **bool**: true or false.
- **vector<pair<int, char>>**: a list where each item is a pair, a number and a letter. For example, [(50,'M'), (2,'I'), (48,'M')]. In code, `.first` gets the number and `.second` gets the letter.
- **const**: "I promise not to change this". The function only reads these inputs.
- **&** (after a type): "use the original, don't make a copy". This is just for speed.
- **\*** (pointer): an **address in memory**, like a locker number. The function builds the form somewhere in memory and gives you the locker number (`bam1_t*`). When you're done, you must empty the locker with `bam_destroy1(...)`. C++ doesn't tidy up for you here; forgetting causes a "memory leak".

### Why two versions with the same name
C++ lets you have two functions with the same name if they take different inputs. This is called **overloading**.

- One version is for **single reads**.
- The other adds extra inputs about the **mate**, for paired reads.

Both call the same helper, `alignment_to_bam_internal`, which does the real work. The single-read version just passes "nothing" values for the mate (empty name, -1, etc.).

## Part 4: What Happens Inside, Step by Step

The code is in `src/alignment.cpp`, in the function `alignment_to_bam_internal`.

### Step 1: Get a blank form
```
bam1_t* bam = bam_init1();
```
This asks the htslib library (the standard BAM toolkit) for an empty record.

### Step 2: Tidy the name
Mates are named `read42/1` and `read42/2`, but the BAM rules say both mates must have the **same** name. So for paired reads, the code chops off the `/1`, `/2`, `_1` or `_2` at the end. Both become `read42`.

### Step 3: Reserve space for the variable-length parts
Some parts of the form have different lengths for different reads: the name, the CIGAR, the letters and the quality scores. BAM stores these four **back to back in one strip of memory**:

`[ name ][ CIGAR ][ letters ][ quality scores ]`

The code first **measures** how much room each part needs, then reserves exactly that much with `calloc`. `calloc` is a function that says "give me this many bytes of memory, all set to zero".

> A **byte** is the basic unit of computer memory. It holds a number from 0 to 255.

### Step 4: Fill in the fixed-size boxes
These boxes are always a single number:

| Box | What gets written | Plain meaning |
|---|---|---|
| pos | refpos | Letter position in the chapter |
| tid | chapter name → number | BAM saves space by storing "chapter #0" instead of "chr1" |
| qual | mapping quality | How confident vg is that this is the right place |
| flag | a set of yes/no answers | Explained just below |
| isize | tlen | For pairs: how long the original DNA piece was |

**Flags: many yes/no answers squeezed into one number.** Think of a row of light switches, where each switch answers one question:

| Switch | Question |
|---|---|
| 1 | Is this read part of a pair? |
| 2 | Does the pair look healthy (same chapter, facing each other, not too far apart)? |
| 4 | Could this read NOT be placed? |
| 8 | Could its mate NOT be placed? |
| 16 | Did it match the other side of the zipper? |
| 32 | Did its mate? |
| 64 / 128 | Is this the first / second of the pair? |
| 256 | Is this a backup placement, not the best one? |

In C++, `flags |= X` means "flip switch X on". The final number is all the "on" switches added together.

### Step 5: Write the name
The name is copied in letter by letter, followed by a few zero bytes. The zeros mark the end of the text and keep things neatly lined up in memory.

### Step 6: Write the CIGAR
Each item, like (50, 'M'), is turned into one compact 4-byte number. A `switch` statement (C++'s "if it's this letter, do this; if that letter, do that") converts each letter into the library's code for it.

If the code sees a letter it doesn't recognise, it **throws an error**: it stops and complains, rather than writing a bad file.

While doing this, it also works out **where the read ends** in the chapter. From that it computes a **bin**, a sort of shelf label that later helps programs quickly find "all reads in this region".

### Step 7: Write the DNA letters, two per byte
There are only a few possible letters, so each one fits in **half a byte**. Two letters share one byte, which saves half the space:

| Letter | A | C | G | T | N (unknown) |
|---|---|---|---|---|---|
| Code | 1 | 2 | 4 | 8 | 15 |

**The zipper-side rule:** BAM always writes letters as they appear on the **front side** of the zipper. If the read matched the **back side** (refrev is true), the code first:

- converts the letters to their **reverse complement** (reverse, then swap A↔T and C↔G);
- **reverses** the quality scores too, so each score still sits next to its letter.

### Step 8: Write the quality scores
These are one byte per letter, copied straight across. If the read has no quality scores, each byte is set to **255**, BAM's way of saying "unknown".

### Step 9: Add optional extras ("tags")
Tags are optional sticky notes on the form, each written as a short name, a type and a value:

| Tag | Meaning |
|---|---|
| AS | The alignment score: how good the fit is (only if the read was placed) |
| RG | The read group: which sample or batch the read came from |
| SS, GR, NR | Extra vg-specific notes, if present |

Any tags the read already had from the original input file are also copied over, converted to the right kind of value (number, text, list…).

### Step 10: Hand back the finished form
```
return bam;
```
Other code then writes it into the output file, and afterwards frees the memory.

## Part 5: One Example, Start to Finish

**vg's notes:**

- name: read42/1 (first of a pair)
- letters: ACCG
- placed on: chr1, position 1000, on the **back side** of the zipper
- match: 4M (all 4 letters line up)
- confidence: 60

**What alignment_to_bam writes:**

| Box | Value | Why |
|---|---|---|
| name | read42 | /1 removed |
| chapter | 0 | chr1 is the first chapter in the table of contents |
| position | 1000 | Copied |
| confidence | 60 | Copied |
| flags | paired + first + back side, … | Switches turned on |
| CIGAR | 4M | Packed into a number |
| letters | CGGT | Reverse complement of ACCG, because it was on the back side |
| quality | reversed | So the scores still match their letters |
| tags | AS:i:… | The score |

## The One-Line Summary

> **alignment_to_bam takes vg's answer about where one piece of DNA belongs and copies it, carefully and in the exact required layout, onto the standard form (BAM) that every other DNA tool can read.**
