/* ----------------------------- MNI Header -----------------------------------
@NAME       : safe_compute_xfm.c
@DESCRIPTION: Routine to calculate a General_transform from a pair of tag
              point lists. This routine is safe in that it returns
              an identity matrix if any error occurs.
@METHOD     : 
@GLOBALS    : 
@CREATED    : April 21, 1994 (Peter Neelin)
@MODIFIED   : $Log: safe_compute_xfm.c,v $
@MODIFIED   : Revision 1.1  1994-11-04 14:45:55  david
@MODIFIED   : Initial revision
@MODIFIED   :
 * Revision 1.1  94/04/22  08:07:20  neelin
 * Initial revision
 * 
---------------------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>
#include <mni.h>
#include <compute_xfm.h>

/* ----------------------------- MNI Header -----------------------------------
@NAME       : safe_compute_transform_from_tags
@INPUT      : npoints - number of pairs of tag points
              tag_list1 - first list of tag points
              tag_list2 - second list of tag points
              trans_type - type of transformation to calculate
@OUTPUT     : transform - computed transform
@RETURNS    : (nothing)
@DESCRIPTION: Routine to calculate a general transform from a pair of lists
              of tag points. The transform is from tag_list2 to tag_list1.
              This routine is safe in that it returns an identity matrix 
              if any error occurs.
@METHOD     : 
@GLOBALS    : 
@CALLS      :
@CREATED    : April 21, 1994 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
public void safe_compute_transform_from_tags(int npoints, 
                                             Real **tag_list1, 
                                             Real **tag_list2, 
                                             Trans_type trans_type,
                                             General_transform *transform)
{
   int fildes[2];
   FILE *fpin, *fpout;
   Status status;
   int statptr;

   /* Create a pipe */
   if (pipe(fildes)) {
      create_linear_transform(transform, NULL);
      return;
   }

   /* Fork */
   if (fork()) {          /* Parent */
      (void) close(fildes[1]);
      fpin = fdopen(fildes[0], "r");
      status = input_transform(fpin, transform);
      (void) fclose(fpin);
      do {
         (void) wait(&statptr);
      } while (WIFSTOPPED(statptr));
      if (WEXITSTATUS(statptr) || status != OK) {
         create_linear_transform(transform, NULL);
         return;
      }
   }

   else {                 /* Child */
      (void) close(fildes[0]);
      fpout = fdopen(fildes[1], "w");
      compute_transform_from_tags(npoints, tag_list1, tag_list2, trans_type,
                                  transform);
      status = output_transform(fpout, NULL, transform);
      (void) fclose(fpout);
      if (status != OK) {
         exit(EXIT_FAILURE);
      }
      else {
         exit(EXIT_SUCCESS);
      }
   }

   return;
}
