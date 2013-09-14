#ifndef DEF_SAM_H
#define DEF_SAM_H

#include "mrf.h"

/**
 *   \file sam.h
 */

// Bitwise flags for FLAGS field in SAM entry
// Notes:
//  1) Flag 0x02, 0x08, 0x20, 0x40, 0x80 are only meaningful when 0x01 is set
//  2) If in a read pair the information on which read is the first in the pair
//     is lost in upstream analysis, flag 0x01 should be present and 0x40 and
//     0x80 are both zero
#define S_READ_PAIRED     0x0001  // Read is paired in sequencing
#define S_PAIR_MAPPED     0x0002  // Read is mapped in proper pair
#define S_QUERY_UNMAPPED  0x0004  // Query sequence is unmapped
#define S_MATE_UNMAPPED   0x0008  // Mate is unmapped
#define S_QUERY_STRAND    0x0010  // Strand of query (0 forward; 1 reverse)
#define S_MATE_STRAND     0x0020  // Strand of mate (0 forward; 1 reverse)
#define S_FIRST           0x0040  // Read is first read in a pair
#define S_SECOND          0x0080  // Read is second read in a pair
#define S_NOT_PRIMARY     0x0100  // Alignment is not primary
#define S_FAILS_CHECKS    0x0200  // Read fails platform/vendor checks
#define S_DUPLICATE       0x0400  // Read is PCR or optical duplicate

typedef enum {
  kCigarAlignmentMatch,   // Alignment match
  kCigarInsertion,        // Insertion into the reference
  kCigarDeletion,         // Deletion from the reference
  kCigarSkipRegion,       // Skipped region from the reference
  kCigarSoftClip,         // Soft clipping
  kCigarHardClip,         // Hard clipping
  kCigarPadding,          // Padding
  kCigarSequenceMatch,    // Sequence match
  kCigarSequenceMismatch, // Sequence mismatch
  kCigarInvalid
} CigarType;

typedef struct {
  CigarType type;
  int length;
} CigarOperation;

/**
 * SamEntry.
 */
typedef struct {
  char *qname;        // Query name
  int flags;          // Bitwise FLAGS field
  char *rname;        // Reference sequence name
  int pos;            // 1-based leftmost position/coordinate of clipped seq.
  int mapq;           // Mapping quality
  char *cigar;        // Extended CIGAR string
  char *mrnm;         // Mate reference sequence name ("=" if same as rname)
  int mpos;           // 1-based leftmost mate position of clipped sequence
  int isize;          // Inferred insert size
  char *seq;          // Query sequence
  char *qual;         // Query quality string
  char *tags;         // Optional tags (actually list, but as string for now)
} SamEntry;

int sortSamEntriesByQname (SamEntry *a, SamEntry *b);
Stringa genCigar (MrfRead *read);
void destroySamEArray (Array a);

void samParser_initFromFile(char* fileName);
void samParser_initFromPipe(char* command);
void samParser_deInit(void);
void samParser_copyEntry(SamEntry **dest, SamEntry *orig);
void samParser_freeEntry(SamEntry *currEntry);
SamEntry* samParser_nextEntry(void);
char* samParser_writeEntry(SamEntry* currSamEntry);
Array samParser_getAllEntries();
Array samParser_getCigar(char* cigar_string);

#endif /* DEF_SAM_H */
