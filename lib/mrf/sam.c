#include <bios/format.h>
#include <bios/linestream.h>
#include <bios/common.h>

#include "sam.h"



/** 
 *   \file sam.c SAM utilities.
 */

static LineStream ls = NULL;


int sortSamEntriesByQname (SamEntry *a, SamEntry *b)
{
  return strcmp (a->qname, b->qname);
}



Stringa genCigar (MrfRead *read)
{
  int ltargetEnd = 0;
  int i;
  Stringa cigar = stringCreate(10);

  for (i = 0; i < arrayMax (read->blocks); i++) {
    MrfBlock *block = arrp (read->blocks, i, MrfBlock);
    if (ltargetEnd > 0)
	  stringAppendf (cigar, "%iN", block->targetStart - ltargetEnd - 1);
    ltargetEnd = block->targetEnd;
    stringAppendf (cigar, "%iM", block->queryEnd - block->queryStart + 1);
  }
  return cigar;
}



void destroySamEArray (Array a)
{
  int i;
  for (i = 0; i < arrayMax(a); i++) {
    SamEntry *currSamE = arrp (a, i, SamEntry);
    free (currSamE->qname);
    free (currSamE->rname);
    free (currSamE->cigar);
    free (currSamE->mrnm);
    if (currSamE->seq)
      free (currSamE->seq);
    if (currSamE->qual)
      free (currSamE->qual);
    if (currSamE->tags)
      free (currSamE->tags);
  }
  arrayDestroy (a);
}




/**
 * Initialize the SAM module from file.
 * @param[in] fileName File name, use "-" to denote stdin
 */
void samParser_initFromFile (char *fileName)
{
  ls = ls_createFromFile (fileName);
  ls_bufferSet (ls,1);
}



/**
 * Initialize the samParser module from pipe.
 * @param[in] command Command to be executed
 */
void samParser_initFromPipe (char *command)
{
  ls = ls_createFromPipe (command);
  ls_bufferSet (ls,1);
}



/**
 * Deinitialize the samParser module.
 */
void samParser_deInit (void)
{
  ls_destroy (ls);
}



/**
 * Deinitialize the SamEntry.
 */
void samParser_freeEntry (SamEntry *currSamEntry) 
{
  if (currSamEntry == NULL) 
    return;
  hlr_free (currSamEntry->qname);
  hlr_free (currSamEntry->rname);
  hlr_free (currSamEntry->cigar);
  hlr_free (currSamEntry->mrnm);
  if (currSamEntry->seq)
    hlr_free (currSamEntry->seq);
  if (currSamEntry->qual)
    hlr_free (currSamEntry->qual);
  if (currSamEntry->tags)
    hlr_free (currSamEntry->tags);
  freeMem (currSamEntry);
}

/**
 * Make a copy of a SamEntry.
 * @pre *dest
 * @post Use samParser_freeEntry to de-allocate the memory
 */
void samParser_copyEntry (SamEntry **dest, SamEntry *orig) 
{
  if( *dest ) freeMem(*dest);
  AllocVar (*dest); 
  (*dest)->qname = hlr_strdup(orig->qname);
  (*dest)->flags = orig->flags;
  (*dest)->rname = hlr_strdup(orig->rname);
  (*dest)->pos   = orig->pos;
  (*dest)->mapq  = orig->mapq;
  (*dest)->cigar = hlr_strdup(orig->cigar);
  (*dest)->mrnm  = hlr_strdup(orig->mrnm);
  (*dest)->mpos  = orig->mpos;
  (*dest)->isize = orig->isize;
  (*dest)->seq   = orig->seq != NULL ? hlr_strdup(orig->seq) : NULL;
  (*dest)->qual  = orig->qual != NULL ? hlr_strdup(orig->qual) : NULL;
  (*dest)->tags  = orig->tags != NULL ? hlr_strdup(orig->tags) : NULL;
}


int isMateUnmapped( SamEntry* samE ) 
{
  if( samE->flags & S_MATE_UNMAPPED )
    return 1;
  else
    return 0;
}

/**
 * Check if it is a valid Sam line.
 * @pre The module has been initialized using samParser_init().
 * @res 0 if not valid, otherwise 1.
 */
int isPaired( SamEntry* samE )
{
  if( samE->flags & S_READ_PAIRED )  
    return 1;
  else 
    return 0;
}

/**
 * Check if it is a valid Sam line.
 * @pre The module has been initialized using samParser_init().
 * @res 0 if not valid, otherwise 1.
 */
int isValidSamLine( Texta tokens ) {
  if (arrayMax (tokens) < 11) {
    textDestroy( tokens );
    return 0;
  } else return 1;
}




static void samParser_processLine (char* line, SamEntry* currSamEntry) 
{
  Texta tokens = NULL;
  int j;
  int hasQual = 0;
  int hasSeqs = 0;
  tokens = textFieldtokP (line, "\t");
  if( isValidSamLine( tokens ) != 1 ) {
    ls_destroy (ls);
    die ("Invalid SAM entry: %s", line);
  }
 
  currSamEntry->qname = hlr_strdup(textItem(tokens, 0));
  currSamEntry->flags = atoi(textItem(tokens, 1));
  currSamEntry->rname = hlr_strdup(textItem(tokens, 2));
  currSamEntry->pos   = atoi(textItem(tokens, 3));
  currSamEntry->mapq  = atoi(textItem(tokens, 4));
  currSamEntry->cigar = hlr_strdup(textItem(tokens, 5));
  currSamEntry->mrnm  = hlr_strdup(textItem(tokens, 6));
  currSamEntry->mpos  = atoi(textItem(tokens, 7));
  currSamEntry->isize = atoi(textItem(tokens, 8));
  currSamEntry->seq   = NULL;
  currSamEntry->qual  = NULL;
  currSamEntry->tags  = NULL;
  // Get tokens
  if (arrayMax (tokens) > 11) {
    Stringa tags = stringCreate (10);
    for (j = 11; j < arrayMax (tokens); j++) {
      if (j > 11)
	stringAppendf (tags, "\t");
      stringAppendf (tags, "%s", textItem (tokens, j));
    }
    currSamEntry->tags = hlr_strdup (string(tags));
    stringDestroy (tags);
  }
  
  if (strcmp (textItem (tokens, 9),  "*") != 0) {
    hasSeqs = 1;
    currSamEntry->seq = hlr_strdup (textItem (tokens, 9));
  }
  if (strcmp (textItem (tokens, 10), "*") != 0) {
    hasQual = 1;
    currSamEntry->qual = hlr_strdup (textItem (tokens, 10));
  } 
  textDestroy( tokens );
}


static SamEntry* samParser_processNextEntry (int freeMemory)
{
  char *line,*pos;
  //char *queryName = NULL;
  //char *prevSamEntryName = NULL;
  static SamEntry *currSamEntry = NULL;

  if (!ls_isEof (ls)) {
    if (freeMemory) {
      samParser_freeEntry (currSamEntry);
      currSamEntry = NULL;
    }
    AllocVar (currSamEntry);
    
    while (line = ls_nextLine (ls)) {
      if (line[0] == '@') {
	continue;
      }
      samParser_processLine (line,currSamEntry); 
      return currSamEntry;
    }
  }
  if (freeMemory) {
    samParser_freeEntry (currSamEntry);
    currSamEntry = NULL;
  }
  return NULL;  
}

/**
 * Read next SAM entry
 * @pre *dest
 * @post Use samParser_freeEntry to de-allocate the memory
 */
SamEntry* samParser_nextEntry (void)
{
  return samParser_processNextEntry (1);
}

/**
 * Writes a SamEntry.
 * @pre The module has been initialized using samParser_init().
 */
char* samParser_writeEntry( SamEntry* currSamEntry) {
  static Stringa buffer = NULL;
  stringCreateClear (buffer,200);
  stringAppendf(buffer, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s", 
		currSamEntry->qname,
		currSamEntry->flags,
		currSamEntry->rname,
		currSamEntry->pos,
		currSamEntry->mapq,
		currSamEntry->cigar,
		currSamEntry->mrnm,
		currSamEntry->mpos,
		currSamEntry->isize,
		currSamEntry->seq,
		currSamEntry->qual,
		currSamEntry->tags);
  return string (buffer);
}
/**
 * Returns an Array of SamEntries.
 * @note The memory belongs to this routine.
 */
Array samParser_getAllEntries ()
{
  Array samQueries;
  SamEntry *currSamEntry;

  samQueries = arrayCreate (100000,SamEntry);
  while (currSamEntry = samParser_processNextEntry (0)) {
    array (samQueries,arrayMax (samQueries),SamEntry) = *currSamEntry;
    freeMem (currSamEntry);
  }
  return samQueries;
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

   Col	Field	Description
1.	QNAME	Query template/pair NAME
2.	FLAG	bitwise FLAG
3.	RNAME	Reference sequence NAME
4.	POS	1-based leftmost POSition/coordinate of clipped sequence
5.	MAPQ	MAPping Quality (Phred-scaled)
6.	CIAGR	extended CIGAR string
7.	MRNM	Mate Reference sequence NaMe (‘=’ if same as RNAME)
8.	MPOS	1-based Mate POSistion
9.	TLEN	inferred Template LENgth (insert size)
10.	SEQ	query SEQuence on the same strand as the reference
11.	QUAL	query QUALity (ASCII-33 gives the Phred base quality)
12+	OPT	variable OPTional fields in the format TAG:VTYPE:VALUE


    Example:

    HWUSI-EAS519_1:5:113:14691:9858	161	chr1	2483	30	54M	=	3349	920	GTGATGCCAGGCATGCCCTTCCCCAGCATCAGGTCTCCAGAGCTGCAGAAGACG	CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCCCCC@CCCCCACACCCCA	MD:Z:54	SM:i:0
HWUSI-EAS519_1:5:113:14691:9858	81	chr1	3349	30	54M	=	2483	-920	GCTGCACCACTGCCTGGCGCTGTGCCCTTCCTTTGCTCTGCCCGCTGGAGACGG	CCCCC@CCCCCCCCCCCBCCCCCCBCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	MD:Z:54	SM:i:3
HWUSI-EAS519_1:5:7:2249:4114	161	chr1	4270	35	54M	=	4593	377	CTGCTCAGTTCTTTATTGATTGGTGTGCCGTTTTCTCTGGAAGCCTCTTAAGAA	CCCCCCCCCCCCCBCCCCCCCCABCBCBCDCCCCCBCCC=CC>CC<C@C@A;CC	MD:Z:54	SM:i:1
HWUSI-EAS519_1:5:18:5639:19185	97	chr1	4270	35	54M	=	4600	384	CTGCTCAGTTCTTTATTGATTGGTGTGCCGGTTTCTCTGGAAGCCTCTTAAGAA	CCCCCCCCCCCCCCCCCCCCCCCDDCCCCC2CCCCCCCCCCCCCC?CCCC@CCC	MD:Z:30T23	SM:i:1
HWUSI-EAS519_1:5:97:7914:12219	97	chr1	4270	36	54M	=	4593	377	CTGCTCAGTTCTTTATTGATTGGTGTGCCGTTTTCTCTGGAAGCCTCTTAAGAA	CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCBCCCCCCCCCCCCCCD	MD:Z:54	SM:i:1
HWUSI-EAS519_1:5:104:8269:13558	97	chr1	4270	33	54M	=	4587	371	CTGCTCAGTTCTTTATTGATTGGTGTGCCGTTTTCTCCGGAAGCCTCTTAAGAA	CCCCCCCCCCCCCCCDCCCCDCCCCCBCCCCCCCC/+(/000000CCCCCC@CC	MD:Z:37T16	SM:i:1

   \endverbatim
 */
