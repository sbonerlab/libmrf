/// @file sam.c
/// @version 0.8.0
/// @since 19 Sep 2013
///
/// @section DESCRIPTION
///
/// Parser for SAM files.

#define _GNU_SOURCE

#include <stdbool.h>

#include <bios/format.h>
#include <bios/log.h>
#include <bios/linestream.h>
#include <bios/common.h>

#include "sam.h"

static LineStream ls = NULL;

int samentry_compare_by_qname(SamEntry *a, SamEntry *b) {
  return strcmp(a->qname, b->qname);
}

bool samentry_is_mate_unmapped(SamEntry* self) {
  return self->flags & S_MATE_UNMAPPED;
}

/**
 * Check if it is a valid Sam line.
 * @pre The module has been initialized using samParser_init().
 * @res 0 if not valid, otherwise 1.
 */
bool samentry_is_paired(SamEntry* self) {
  return self->flags & S_READ_PAIRED;
}

/**
 * Make a copy of a SamEntry.
 * @pre *dest
 * @post Use samParser_freeEntry to de-allocate the memory
 */
void samentry_copy(SamEntry **dest, SamEntry *orig) {
  if (*dest) {
    free(*dest);
  }
  *dest = samentry_new();
  (*dest)->qname = strdup(orig->qname);
  (*dest)->flags = orig->flags;
  (*dest)->rname = strdup(orig->rname);
  (*dest)->pos   = orig->pos;
  (*dest)->mapq  = orig->mapq;
  (*dest)->cigar = strdup(orig->cigar);
  (*dest)->mrnm  = strdup(orig->mrnm);
  (*dest)->mpos  = orig->mpos;
  (*dest)->isize = orig->isize;
  (*dest)->seq   = (orig->seq != NULL) ? strdup(orig->seq) : NULL;
  (*dest)->qual  = (orig->qual != NULL) ? strdup(orig->qual) : NULL;
  (*dest)->tags  = (orig->tags != NULL) ? strdup(orig->tags) : NULL;
}

/**
 * Writes a SamEntry.
 * @pre The module has been initialized using samParser_init().
 */
char* samentry_to_string(SamEntry* self) {
  static Stringa buffer = NULL;
  stringCreateClear(buffer,200);
  stringAppendf(buffer, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s", 
                self->qname,
                self->flags,
                self->rname,
                self->pos,
                self->mapq,
                self->cigar,
                self->mrnm,
                self->mpos,
                self->isize,
                self->seq,
                self->qual,
                self->tags);
  return string(buffer);
}

bool samentry_same_name(SamEntry* query, SamEntry* mate, char delim) {
  char* pos = strchr(query->qname, delim);
  if (*pos != '\0') {
    *pos='\0';
  }
  pos = strchr(mate->qname, delim);
  if (*pos != '\0') {
    *pos='\0';
  }
  return strcmp(query->qname, mate->qname) == 0;
}

SamEntry* samentry_new(void) {
  SamEntry* entry = (SamEntry*) malloc(sizeof(*entry));
  if (entry == NULL) {
    return NULL;
  }
  memset(entry, 0, sizeof(entry));
  return entry;
}

/**
 * Free the SamEntry.
 */
void samentry_free(SamEntry* self) {
  if (self == NULL) {
    return;
  }
  free(self->qname);
  free(self->rname);
  free(self->cigar);
  if (self->cigar_ops != NULL) {
    arrayDestroy(self->cigar_ops);
  }
  free(self->mrnm);
  if (self->seq != NULL) {
    free(self->seq);
  }
  if (self->qual != NULL) {
    free(self->qual);
  }
  if (self->tags != NULL) {
    free(self->tags);
  }
  free(self);
}

static inline CigarType samparser_get_cigar_type(char c) {
  switch (c) {
  case 'M': 
    return kCigarAlignmentMatch;
  case 'I': 
    return kCigarInsertion;
  case 'D': 
    return kCigarDeletion;
  case 'N': 
    return kCigarSkippedRegion;
  case 'S': 
    return kCigarSoftClipping;
  case 'H':
    return kCigarHardClipping;
  case 'P':
    return kCigarPadding;
  case '=':
    return kCigarSequenceMatch;
  case 'X':
    return kCigarSequenceMismatch;
  default:
    return kCigarInvalid;
  }
}

Array samparser_parse_cigar(char* cigar_string) {
  Array cigar_operations = arrayCreate(5, CigarOperation);
  Texta tokens = textFieldtokP(cigar_string, "MNSDIP");
  size_t j = 0;
  for (int i = 0; i < arrayMax(tokens) - 1; ++i) {
    while (isdigit(cigar_string[j]) && cigar_string[j] != '\0') {
      ++j;
    }
    CigarType type = samparser_get_cigar_type(cigar_string[j]);
    CigarOperation operation;
    operation.type = type;
    operation.length = atoi(textItem(tokens, i));
    array(cigar_operations, arrayMax(cigar_operations), CigarOperation) =
        operation;
    if (cigar_string[j] != '\0') {
      ++j;
    }
  }
  return cigar_operations;
}

Stringa samparser_mrfread_to_cigar(MrfRead* read) {
  int ltargetEnd = 0;
  Stringa cigar = stringCreate(10);
  for (int i = 0; i < arrayMax(read->blocks); i++) {
    MrfBlock *block = arrp(read->blocks, i, MrfBlock);
    if (ltargetEnd > 0) {
      stringAppendf(cigar, "%iN", block->targetStart - ltargetEnd - 1);
    }
    ltargetEnd = block->targetEnd;
    stringAppendf(cigar, "%iM", block->queryEnd - block->queryStart + 1);
  }
  return cigar;
}

/**
 * Initialize the SAM module from file.
 * @param[in] fileName File name, use "-" to denote stdin
 */
SamParser* samparser_from_file(char* filename) {
  SamParser* parser = (SamParser*) malloc(sizeof(*parser));
  if (parser == NULL) {
    return NULL;
  }
  parser->ls = ls_createFromFile(filename);
  ls_bufferSet(parser->ls, 1);
  return parser;
}

/**
 * Initialize the samParser module from pipe.
 * @param[in] command Command to be executed
 */
SamParser* samparser_from_pipe(char* command) {
  SamParser* parser = (SamParser*) malloc(sizeof(*parser));
  if (parser == NULL) {
    return NULL;
  }
  parser->ls = ls_createFromPipe(command);
  ls_bufferSet(parser->ls, 1);
  return parser;
}

/**
 * Deinitialize the samParser module.
 */
void samparser_free(SamParser* self) {
  ls_destroy(self->ls);
  free(self);
}

/**
 * Check if it is a valid Sam line.
 * @pre The module has been initialized using samParser_init().
 * @res 0 if not valid, otherwise 1.
 */
bool samparser_valid_line(Texta tokens) {
  if (arrayMax(tokens) < 11) {
    textDestroy(tokens);
    return false;
  } else {
    return true;
  }
}

static void samparser_process_line(SamParser* self, char* line, 
                                   SamEntry* entry) {
  int hasQual = 0;
  int hasSeqs = 0;
  Texta tokens = textFieldtokP(line, "\t");
  if (samparser_valid_line(tokens) == false) {
    ls_destroy(self->ls);
    die("Invalid SAM entry: %s", line);
  }
 
  entry->qname = strdup(textItem(tokens, 0));
  entry->flags = atoi(textItem(tokens, 1));
  entry->rname = strdup(textItem(tokens, 2));
  entry->pos   = atoi(textItem(tokens, 3));
  entry->mapq  = atoi(textItem(tokens, 4));
  entry->cigar = strdup(textItem(tokens, 5));
  entry->mrnm  = strdup(textItem(tokens, 6));
  entry->mpos  = atoi(textItem(tokens, 7));
  entry->isize = atoi(textItem(tokens, 8));
  entry->seq   = NULL;
  entry->qual  = NULL;
  entry->tags  = NULL;

  // Parse seq, qual, and tags from remaining tokens.
  if (arrayMax(tokens) > 11) {
    Stringa tags = stringCreate(10);
    for (int j = 11; j < arrayMax(tokens); j++) {
      if (j > 11) {
        stringAppendf(tags, "\t");
      }
      stringAppendf(tags, "%s", textItem(tokens, j));
    }
    entry->tags = hlr_strdup(string(tags));
    stringDestroy(tags);
  }
  
  if (strcmp(textItem(tokens, 9),  "*") != 0) {
    hasSeqs = 1;
    entry->seq = hlr_strdup(textItem(tokens, 9));
  }
  if (strcmp(textItem(tokens, 10), "*") != 0) {
    hasQual = 1;
    entry->qual = hlr_strdup(textItem(tokens, 10));
  } 
  textDestroy(tokens);

  // Parse cigar
  entry->cigar_ops = samparser_parse_cigar(entry->cigar);
}

/**
 * Returns a pointer to next SamEntry. 
 * @pre The module has been initialized using samParser_init().
 * Parse entries of the following format:
   \verbatim
  
   Output (obtained from running sam -h)
   ----------------------------------------

   The 'sam' aligner outputs each alignment on a separate line.  Each
   line is a collection of 8 fields separated by tabs; from left to
   right, the fields are:

   Col   Field   Description
   1.    QNAME   Query template/pair NAME
   2.    FLAG    bitwise FLAG
   3.    RNAME   Reference sequence NAME
   4.    POS     1-based leftmost POSition/coordinate of clipped sequence
   5.    MAPQ    MAPping Quality (Phred-scaled)
   6.    CIAGR   extended CIGAR string
   7.    MRNM    Mate Reference sequence NaMe (‘=’ if same as RNAME)
   8.    MPOS    1-based Mate POSistion
   9.    TLEN    inferred Template LENgth (insert size)
   10.   SEQ     query SEQuence on the same strand as the reference
   11.   QUAL    query QUALity (ASCII-33 gives the Phred base quality)
   12.   OPT     variable OPTional fields in the format TAG:VTYPE:VALUE

   Example:

   HWUSI-EAS519_1:5:113:14691:9858  161  chr1  2483  30  54M  =  3349  920  GTGATGCCAGGCATGCCCTTCCCCAGCATCAGGTCTCCAGAGCTGCAGAAGACG  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCC@CCCCCACACCCCA  MD:Z:54     SM:i:0
   HWUSI-EAS519_1:5:113:14691:9858  81   chr1  3349  30  54M  =  2483  -920 GCTGCACCACTGCCTGGCGCTGTGCCCTTCCTTTGCTCTGCCCGCTGGAGACGG  CCCCC@CCCCCCCCCCCBCCCCCCBCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  MD:Z:54     SM:i:3
   HWUSI-EAS519_1:5:7:2249:4114     161  chr1  4270  35  54M  =  4593  377  CTGCTCAGTTCTTTATTGATTGGTGTGCCGTTTTCTCTGGAAGCCTCTTAAGAA  CCCCCCCCCCCCCBCCCCCCCCABCBCBCDCCCCCBCCC=CC>CC<C@C@A;CC  MD:Z:54     SM:i:1
   HWUSI-EAS519_1:5:18:5639:19185   97   chr1  4270  35  54M  =  4600  384  CTGCTCAGTTCTTTATTGATTGGTGTGCCGGTTTCTCTGGAAGCCTCTTAAGAA  CCCCCCCCCCCCCCCCCCCCCCCDDCCCCC2CCCCCCCCCCCCCC?CCCC@CCC  MD:Z:30T23  SM:i:1
   HWUSI-EAS519_1:5:97:7914:12219   97   chr1  4270  36  54M  =  4593  377  CTGCTCAGTTCTTTATTGATTGGTGTGCCGTTTTCTCTGGAAGCCTCTTAAGAA  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCBCCCCCCCCCCCCCCD  MD:Z:54     SM:i:1
   HWUSI-EAS519_1:5:104:8269:13558  97   chr1  4270  33  54M  =  4587  371  CTGCTCAGTTCTTTATTGATTGGTGTGCCGTTTTCTCCGGAAGCCTCTTAAGAA  CCCCCCCCCCCCCCCDCCCCDCCCCCBCCCCCCCC/+(/000000CCCCCC@CC  MD:Z:37T16  SM:i:1

   \endverbatim
 */
static SamEntry* samparser_process_next_entry(SamParser* self, bool free_memory) {
  static SamEntry *current_entry = NULL;
  if (!ls_isEof (ls)) {
    if (free_memory) {
      samentry_free(current_entry);
      current_entry = NULL;
    }
    current_entry = samentry_new();
    char* line = NULL;
    while (line = ls_nextLine(self->ls)) {
      if (line[0] == '@') {
        continue;
      }
      samparser_process_line(self, line, current_entry); 
      return current_entry;
    }
  }
  if (free_memory != false) {
    samentry_free(current_entry);
    current_entry = NULL;
  }
  return NULL;  
}

/**
 * Read next SAM entry
 * @pre *dest
 * @post Use samParser_freeEntry to de-allocate the memory
 */
SamEntry* samparser_next_entry(SamParser* self) {
  return samparser_process_next_entry(self, 0);
}

/**
 * Returns an Array of SamEntries.
 * @note The memory belongs to this routine.
 */
Array samparser_get_all_entries(SamParser* parser) {
  Array samQueries = arrayCreate (100000,SamEntry);
  SamEntry* entry;
  while (entry = samparser_process_next_entry(parser, 0)) {
    array(samQueries, arrayMax(samQueries), SamEntry) = *entry;
    freeMem(entry);
  }
  return samQueries;
}
