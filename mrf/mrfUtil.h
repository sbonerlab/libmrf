/// @file mrfUtil.c 
/// @version 0.8.0
/// @since 19 Aug 2013
///
/// @section DESCRIPTION
///
/// MRF utilities

#ifndef DEF_MRF_UTIL_H
#define DEF_MRF_UTIL_H

/// @brief Structure representing a TAR.
typedef struct {
  char* targetName;
  int start;
  int end;
} Tar;

extern Array readTarsFromBedFile(char *fileName);

#endif /* DEF_MRF_UTIL_H */
