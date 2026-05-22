
/******************************************************************************
 *
 *  This file is part of meryl-utility, a collection of miscellaneous code
 *  used by Meryl, Canu and others.
 *
 *  This software is based on:
 *    'Canu' v2.0              (https://github.com/marbl/canu)
 *  which is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "sequence.H"
#include "files.H"

#include "htslib/hts/sam.h"

#include <vector>
#include <string>

using namespace merylutil::sequence;

//
//  A tiny "expected record" container we build alongside the file we write
//  so the verification step can compare byte-for-byte.
//
struct expRec {
  std::string  ident;
  std::string  bases;
  std::string  quals;   //  Empty for FASTA, encoded ASCII '!'..'~' for FASTQ.
};

static uint32  g_fails = 0;

#define CHECK(cond, msg)                                                  \
  do {                                                                    \
    if (!(cond)) {                                                        \
      fprintf(stderr, "FAIL: %s  (%s:%d)\n", msg, __FILE__, __LINE__);    \
      g_fails++;                                                          \
    }                                                                     \
  } while (0)


//
//  Build a deterministic ACGT sequence of length 'len' (and matching qvs).
//
static
void
makeSeq(std::string &s, std::string &q, uint32 len, uint32 salt) {
  static const char  bases[4] = { 'A', 'C', 'G', 'T' };
  s.resize(len);
  q.resize(len);
  for (uint32 i=0; i<len; i++) {
    s[i] = bases[(i + salt) & 0x3];
    q[i] = '!' + ((i * 7 + salt) % 40);   //  0..39 -> '!'..'H'
  }
}


//
//  Write a FASTA file containing the records in 'recs'.
//  'wrap' is the line-break length; 0 means single-line (no wrapping).
//
static
void
writeFASTA(char const *filename, std::vector<expRec> &recs, uint64 wrap) {
  FILE *F = merylutil::openOutputFile(filename);
  for (size_t i=0; i<recs.size(); i++)
    merylutil::outputFASTA(F,
                           recs[i].bases.c_str(),
                           recs[i].bases.size(),
                           wrap,
                           "%s", recs[i].ident.c_str());
  merylutil::closeFile(F, filename);
}


//
//  Write a FASTQ file containing the records in 'recs'.
//  Qualities are passed as the ASCII-encoded form ('char *' overload).
//
static
void
writeFASTQ(char const *filename, std::vector<expRec> &recs) {
  FILE *F = merylutil::openOutputFile(filename);
  for (size_t i=0; i<recs.size(); i++)
    merylutil::outputFASTQ(F,
                           recs[i].bases.c_str(),
                           recs[i].quals.c_str(),
                           recs[i].bases.size(),
                           "%s", recs[i].ident.c_str());
  merylutil::closeFile(F, filename);
}


//
//  Read 'filename' with loadSequence() and confirm every record matches recs.
//
static
void
verify_loadSequence(char const          *filename,
                    std::vector<expRec> &recs,
                    bool                 hasQuals,
                    char const          *tag) {
  dnaSeqFile *F = openSequenceFile(filename);
  CHECK(F != nullptr, "openSequenceFile failed in verify_loadSequence");
  if (F == nullptr)  return;

  dnaSeq  S;
  size_t  n = 0;

  while (F->loadSequence(S)) {
    char msg[256];

    if (n >= recs.size()) {
      snprintf(msg, sizeof(msg), "[%s] loadSequence returned more records than expected (%zu)", tag, n+1);
      CHECK(false, msg);
      n++;
      continue;
    }

    snprintf(msg, sizeof(msg), "[%s] record %zu ident mismatch ('%s' vs '%s')",
             tag, n, S.ident(), recs[n].ident.c_str());
    CHECK(strcmp(S.ident(), recs[n].ident.c_str()) == 0, msg);

    snprintf(msg, sizeof(msg), "[%s] record %zu length mismatch (%lu vs %zu)",
             tag, n, (unsigned long)S.length(), recs[n].bases.size());
    CHECK(S.length() == recs[n].bases.size(), msg);

    snprintf(msg, sizeof(msg), "[%s] record %zu bases mismatch", tag, n);
    CHECK(strcmp(S.bases(), recs[n].bases.c_str()) == 0, msg);

    if (hasQuals) {
      bool qok = true;
      for (uint64 i=0; i<S.length(); i++) {
        //  Stored quals are integer (0..93); recs[].quals is ASCII ('!'+v).
        if ((uint8)(S.quals()[i] + '!') != (uint8)recs[n].quals[i]) {
          qok = false;
          break;
        }
      }
      snprintf(msg, sizeof(msg), "[%s] record %zu quals mismatch", tag, n);
      CHECK(qok, msg);
    }
    n++;
  }

  {
    char msg[128];
    snprintf(msg, sizeof(msg), "[%s] loadSequence returned %zu records, expected %zu",
             tag, n, recs.size());
    CHECK(n == recs.size(), msg);
  }

  delete F;
}


//
//  Read 'filename' with chunked loadBases() and verify the concatenated bases
//  for each record match recs.  bufSize is the per-call maxLength.
//
static
void
verify_loadBases(char const          *filename,
                 std::vector<expRec> &recs,
                 uint64               bufSize,
                 char const          *tag) {
  dnaSeqFile *F = openSequenceFile(filename);
  CHECK(F != nullptr, "openSequenceFile failed in verify_loadBases");
  if (F == nullptr)  return;

  char   *buf = new char [bufSize];
  size_t  n   = 0;

  std::string  acc;

  while (true) {
    uint64  seqLen = 0;
    bool    eos    = false;

    bool ok = F->loadBases(buf, bufSize, seqLen, eos);

    //  EOF with nothing buffered: we're done.  Don't count this as a record;
    //  loadBases() leaves endOfSequence==true at EOF even though no record
    //  was returned.
    if (!ok && seqLen == 0 && acc.empty())
      break;
 
    if (seqLen > 0)
      acc.append(buf, buf + seqLen);

    if (eos) {
      char msg[256];

      if (n >= recs.size()) {
        snprintf(msg, sizeof(msg), "[%s,buf=%lu] loadBases returned more records than expected (%zu)",
                 tag, (unsigned long)bufSize, n+1);
        CHECK(false, msg);
      }
      else {
        snprintf(msg, sizeof(msg),
                 "[%s,buf=%lu] record %zu loadBases bases mismatch (got len %zu, expected len %zu)",
                 tag, (unsigned long)bufSize, n, acc.size(), recs[n].bases.size());
        CHECK(acc == recs[n].bases, msg);
      }
      n++;
      acc.clear();
    }

    if (!ok)
      break;
  }

  {
    char msg[160];
    snprintf(msg, sizeof(msg),
             "[%s,buf=%lu] loadBases produced %zu records, expected %zu",
             tag, (unsigned long)bufSize, n, recs.size());
    CHECK(n == recs.size(), msg);
  }

  delete [] buf;
  delete F;
}


//
//  Build a small BAM file with 'recs' (forward strand only so the reader
//  doesn't reverse-complement them).
//
static
void
writeBAM(char const *filename, std::vector<expRec> &recs) {
  htsFile *w = hts_open(filename, "wb");
  assert(w != nullptr);

  sam_hdr_t *hdr = bam_hdr_init();
  assert(hdr != nullptr);
  int r = sam_hdr_add_line(hdr, "SQ", "SN", "t1", "LN", "100000", NULL);
  assert(r == 0);
  r = sam_hdr_write(w, hdr);
  assert(r == 0);

  bam1_t *b = bam_init1();
  for (size_t i=0; i<recs.size(); i++) {
    //  Unmapped, forward strand (no BAM_FREVERSE -> no rev-comp on read).
    r = bam_set1(b,
                 recs[i].ident.size(), recs[i].ident.c_str(),
                 BAM_FUNMAP, -1, -1, 0,
                 0, NULL,                       //  no cigar
                 -1, -1, 0,                     //  no mate
                 recs[i].bases.size(),
                 recs[i].bases.c_str(),
                 recs[i].quals.c_str(),         //  ASCII '!'+v
                 0);
    assert(r >= 0);
    r = sam_write1(w, hdr, b);
    assert(r >= 0);
  }
  bam_destroy1(b);
  sam_hdr_destroy(hdr);
  r = hts_close(w);
  assert(r == 0);
}


//
//  For BAM verification, build the expected ASCII-uppercase bases (the
//  hts reader produces only [ACGTN] uppercase via its decode table).
//
static
void
expectedBAMbases(std::vector<expRec> &in, std::vector<expRec> &out) {
  out = in;
  for (size_t i=0; i<out.size(); i++)
    for (size_t j=0; j<out[i].bases.size(); j++) {
      char c = out[i].bases[j];
      if      (c == 'a') c = 'A';
      else if (c == 'c') c = 'C';
      else if (c == 'g') c = 'G';
      else if (c == 't') c = 'T';
      else if (c == 'n') c = 'N';
      out[i].bases[j] = c;
    }
}


int
main(int argc, char **argv) {

  //
  //  Build a set of records:
  //   - single-line (length < wrap),
  //   - exact multiple of wrap,
  //   - not a multiple of wrap (this is the tricky multi-line case).
  //
  std::vector<expRec>  recs(3);

  recs[0].ident = "r0_short";    makeSeq(recs[0].bases, recs[0].quals,  7, 0);
  recs[1].ident = "r1_exact";    makeSeq(recs[1].bases, recs[1].quals, 30, 1);
  recs[2].ident = "r2_offcut";   makeSeq(recs[2].bases, recs[2].quals, 53, 2);

  uint64 const  wrap = 10;

  //
  //  Multi-line FASTA round-trip.
  //
  {
    char const *fn = "sequenceReadTest.ml.fasta";
    writeFASTA(fn, recs, wrap);
    verify_loadSequence(fn, recs, /*hasQuals*/false, "FASTA-multiline");
    verify_loadBases   (fn, recs, /*buf*/    7, "FASTA-multiline");
    verify_loadBases   (fn, recs, /*buf*/   16, "FASTA-multiline");
    verify_loadBases   (fn, recs, /*buf*/65536, "FASTA-multiline");
    if (g_fails == 0)  merylutil::unlink(fn);
  }

  //
  //  Single-line FASTA round-trip.
  //
  {
    char const *fn = "sequenceReadTest.sl.fasta";
    writeFASTA(fn, recs, /*wrap*/0);
    verify_loadSequence(fn, recs, /*hasQuals*/false, "FASTA-singleline");
    verify_loadBases   (fn, recs, /*buf*/    7, "FASTA-singleline");
    verify_loadBases   (fn, recs, /*buf*/65536, "FASTA-singleline");
    if (g_fails == 0)  merylutil::unlink(fn);
  }

  //
  //  FASTQ round-trip.
  //
  {
    char const *fn = "sequenceReadTest.fastq";
    writeFASTQ(fn, recs);
    verify_loadSequence(fn, recs, /*hasQuals*/true, "FASTQ");
    verify_loadBases   (fn, recs, /*buf*/    7, "FASTQ");
    verify_loadBases   (fn, recs, /*buf*/65536, "FASTQ");
    if (g_fails == 0)  merylutil::unlink(fn);
  }

  //
  //  BAM round-trip.
  //
  {
    char const *fn = "sequenceReadTest.bam";
    writeBAM(fn, recs);

    std::vector<expRec> expBAM;
    expectedBAMbases(recs, expBAM);

    verify_loadSequence(fn, expBAM, /*hasQuals*/false, "BAM");
    verify_loadBases   (fn, expBAM, /*buf*/    7, "BAM");
    verify_loadBases   (fn, expBAM, /*buf*/65536, "BAM");
    if (g_fails == 0)  merylutil::unlink(fn);
  }

  if (g_fails == 0) {
    fprintf(stderr, "Success!\n");
    return 0;
  }

  fprintf(stderr, "FAILED with %u check(s).\n", g_fails);
  return 1;
}
