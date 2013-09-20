/// @file segmentationUtil.h
/// @version 0.8.0
/// @since 19 Aug 2013
///
/// @section DESCRIPTION
///
/// Segmentation utilities.

#ifndef DEF_SEGMENTATION_UTIL_H
#define DEF_SEGMENTATION_UTIL_H

typedef struct {
  char* targetName;
  int start;
  int end;
} Tar;

typedef struct {
  int position;
  float value;
} Wig;

void performSegmentation(Array tars, Array wigs, char* targetName, 
                         double threshold, int maxGap, int minRun);

#endif
