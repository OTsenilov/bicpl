#ifndef  _VOLUME_STUFF_H
#define  _VOLUME_STUFF_H

/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#ifndef lint
static char vols_rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Include/bicpl/vols.h,v 1.1 2000-02-06 15:30:37 stever Exp $";
#endif

#include  <volume_io.h>
#include  <bicpl/colour_coding.h>

typedef  enum  { FOUR_NEIGHBOURS, EIGHT_NEIGHBOURS } Neighbour_types;

typedef struct
{
    int                    x, y;
    Volume                 src_volume;
    Volume                 dest_volume;
    General_transform      transform;
} resample_struct;

#ifndef  public
#define       public   extern
#define       public_was_defined_here
#endif

#include  <bicpl/vol_prototypes.h>

#ifdef  public_was_defined_here
#undef       public
#undef       public_was_defined_here
#endif

#endif