/// @file mrfUtil.c 
/// @version 0.8.0
/// @since 19 Aug 2013
///
/// @section DESCRIPTION
///
/// MRF utilities

#include <bios/format.h>
#include <bios/log.h>
#include <bios/linestream.h>

#include "mrfUtil.h"

Array readTarsFromBedFile (char *fileName) {
  Array tars;
  Tar *currTar;
  LineStream ls;
  WordIter w;
  char *line;
 
  tars = arrayCreate (100000,Tar);
  ls = ls_createFromFile (fileName);
  while (line = ls_nextLine (ls)) {
    if (strStartsWithC (line,"browser") || strStartsWithC (line,"track")) {
      continue;
    }
    w = wordIterCreate (line,"\t",0);
    currTar = arrayp (tars,arrayMax (tars),Tar);
    currTar->targetName = hlr_strdup (wordNext (w));
    currTar->start = atoi (wordNext (w));
    currTar->end = atoi (wordNext (w));
    wordIterDestroy (w);
  }
  ls_destroy (ls);
  return tars;
}
