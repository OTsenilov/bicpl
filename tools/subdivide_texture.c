/*
   SUBDIVIDE_TEXTURE

   Extend a .txt file to the next refinement level. Use nearest neighbour.

   subdivide_texture 

   Values: in.obj = input object file on fine level
           in.txt = input texture on coarse level
           out.txt = output texture on fine level

   By: Claude Lepage, October 2010.

   COPYRIGHT: McConnell Brain Imaging Center, 
              Montreal Neurological Institute,
              Department of Psychology,
              McGill University, Montreal, Quebec, Canada. 
  
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
*/

#include <math.h>
#include <stdio.h>
#include <volume_io.h>
#include <bicpl.h>

// Prototypes of functions in this file.

static void usage( char * );
static Status read_surface_obj( STRING, int *, Point *[],
                                Vector *[], int *, int *[], int *[], int **[] );
static Status get_surface_neighbours( polygons_struct *, int *[],
                                      int ** [] );
static Real * read_scalar( int n_points, char * in_file );
static void save_scalar( int n_points, Real scalar[], char * out_file );


// Main program.

int main( int argc, char * argv[] ) {

  int      i, j, k, jj, kk, pp, v1, v2, opp1, opp2, t1, t2, count, found3;
  int      target_nodes, changed, num_swapped, histo[101];
  Real     thresholdlen, minlen, maxlen, factor;

  int      n_points;           // number of grid points per object
  int      old_n_points;       // number of grid points per object
  int      n_elems;            // number of triangles per object
  Point  * coords;             // coordinates
  Vector * normals;            // normal vectors
  int    * connec;             // connectivity
  int    * n_ngh = NULL;       // node neighbours (inverse connectivity)
  int   ** ngh = NULL;

  FILE   * fp;

  if( argc != 4 ) {
    usage( argv[0] );
    return( 1 );
  }

  // Read the surface file.
  if( read_surface_obj( argv[1], &n_points, &coords, &normals,
                        &n_elems, &connec, &n_ngh, &ngh ) != OK ) {
    return 1;
  }
  FREE( coords );
  FREE( normals );
  FREE( connec );

  // Read the texture file.
  old_n_points = n_points - ( ( n_elems / 4 ) * 3 ) / 2;
  Real * labels = read_scalar( old_n_points, argv[2] );

  // Process the texture map: interpolate labels to new labels.

  Real * newlabels = (Real *)malloc( n_points * sizeof( Real ) );
  for( i = 0; i < old_n_points; i++ ) {
    newlabels[i] = labels[i];
  }

  for( i = old_n_points; i < n_points; i++ ) {
    Real val = -1;
    int found = 0;
    for( k = 0; k < n_ngh[i]; k++ ) {
      int p1 = ngh[i][k];
      if( p1 < old_n_points ) {
        if( !found ) {
          val = labels[p1];
        } else {
          if( found == 1 ) {
            if( val < labels[p1] ) {   // take the minimum label.
              newlabels[i] = val;
            } else {
              newlabels[i] = labels[p1];
            }
          } else {
            printf( "Too many old neighbours at node %d\n", i );
            exit(1);
          }
        }
        found++;
      }
    }
  }

  // Save the new texture file.
  save_scalar( n_points, newlabels, argv[3] );
  free( labels );
  free( newlabels );

  return 0;
}

// -------------------------------------------------------------------
// Help message on how to use this module.
//
static void usage( char * executable_name ) {

  STRING  usage_format = "\
Usage: %s in.obj data.txt out.obj\n\
Values: in.obj = input object file\n\
        out.obj = output object file\n\
        sphere.obj = stereotaxic sphere model\n\n";

  print_error( usage_format, executable_name );
}


// -------------------------------------------------------------------
// Load the cortical surface.
//
// filename: name of the .obj file
// n_points: the number of the vertices
// points: (x,y,z) coordinates
// normals: normal vectors
// n_elem: number of triangles
// connec: connectivity of triangles
// n_neighbours: number of vertices around each node
// neighbours: the set of ordered triangle consisting of the vertices
//
static Status read_surface_obj( STRING filename,
                                 int * n_points,
                                 Point * points[],
                                 Vector * normals[],
                                 int * n_elem,
                                 int * connec[],
                                 int * n_neighbours[],
                                 int ** neighbours[] ) {

  int               i, n_objects;
  object_struct  ** object_list;
  polygons_struct * surface;
  File_formats      format;
  STRING            expanded;

  expanded = expand_filename( filename );   // why?????

  int err = input_graphics_file( expanded, &format, &n_objects,
                                 &object_list );

  if( err != OK ) {
    print_error( "Error reading file %s\n", expanded );
    return( ERROR );
  }

  if( n_objects != 1 || 
      ( n_objects == 1 && get_object_type(object_list[0]) != POLYGONS ) ) {
    print_error( "Error in contents of file %s\n", expanded );
    return( ERROR );
  }

  delete_string( expanded );

  surface = get_polygons_ptr( object_list[0] );

  int ntri = 0, nquad = 0, unknown = 0;
  int start_ind = 0;
  for( i = 0; i < surface->n_items; i++ ) {
    int nn = surface->end_indices[i] - start_ind;
    start_ind = surface->end_indices[i];
    if( nn == 3 ) {
      ntri++;
    } else {
     if( nn == 4 ) {
       nquad++;
     } else {
       unknown++;
       printf( "face with %d nodes\n", nn );
     }
   }
  }
  printf( "%d triangles, %d quads, %d unknown faces in mesh\n", ntri, nquad, unknown );

  // Check if all polygons are triangles.

  if( 3 * surface->n_items != surface->end_indices[surface->n_items-1] ) {
    printf( "Error: Surface must contain only triangular polygons.\n" );
    delete_object_list( n_objects, object_list );
    return ERROR;
  }

  // Make a copy of the coordinates, the normals, and the
  // connectivity since delete_object_list will destroy them.

  *n_points = surface->n_points;
  *n_elem = surface->n_items;
  ALLOC( *points, surface->n_points );
  ALLOC( *normals, surface->n_points );

  for( i = 0; i < *n_points; i++ ) {
    (*points)[i].coords[0] = surface->points[i].coords[0];
    (*points)[i].coords[1] = surface->points[i].coords[1];
    (*points)[i].coords[2] = surface->points[i].coords[2];
    (*normals)[i].coords[0] = surface->normals[i].coords[0];
    (*normals)[i].coords[1] = surface->normals[i].coords[1];
    (*normals)[i].coords[2] = surface->normals[i].coords[2];
  }

  if( connec ) {
    ALLOC( *connec, 3*surface->n_items );
    for( i = 0; i < 3*surface->n_items; i++ ) {
      (*connec)[i] = surface->indices[i];
    }
  }

  if( n_neighbours && neighbours ) {
    get_surface_neighbours( surface, n_neighbours, neighbours );
  }

  delete_object_list( n_objects, object_list );

  return( OK );
}


// -------------------------------------------------------------------
// Construct the edges around each node. The edges are sorted to
// make an ordered closed loop.
//
private Status get_surface_neighbours( polygons_struct * surface,
                                       int * n_neighbours_return[],
                                       int ** neighbours_return[] ) {

  int    i, j, k, jj;
  int  * tri;
  int  * n_ngh;
  int ** ngh;
  int  * ngh_array;

  // Check if all polygons are triangles.

  if( 3 * surface->n_items != surface->end_indices[surface->n_items-1] ) {
    printf( "Surface must contain only triangular polygons.\n" );
    return ERROR;
  }

  // Check if the node numbering starts at 0 or 1.

  int min_idx, max_idx;

  min_idx = 100*surface->n_points;  // anything big
  max_idx = 0;                      // anything small

  for( i = 0; i < 3*surface->n_items; i++ ) {
    if( surface->indices[i] < min_idx ) min_idx = surface->indices[i];
    if( surface->indices[i] > max_idx ) max_idx = surface->indices[i];
  }

  // Shift numbering to start at zero, for array indexing. Note
  // that we don't care if surface->indices array is modified.

  if( min_idx != 0 ) {
    for( i = 0; i < 3*surface->n_items; i++ ) {
      surface->indices[i] -= min_idx;
    }
  }

  // Count number of triangles attached to each node.

  ALLOC( n_ngh, surface->n_points );
  ALLOC( ngh, surface->n_points );
  ALLOC( ngh_array, 3*surface->n_items );

  for( i = 0; i < surface->n_points; i++ ) {
    n_ngh[i] = 0;
  }

  for( i = 0; i < 3*surface->n_items; i++ ) {
    n_ngh[surface->indices[i]]++;
    ngh_array[i] = -1;
  }

  int max_ngh = 0;
  int sum_ngh = 0;
  for( i = 0; i < surface->n_points; i++ ) {
    ngh[i] = &(ngh_array[sum_ngh]);
    sum_ngh += n_ngh[i];
    max_ngh = MAX( max_ngh, n_ngh[i] );
  }

  // At first, store the indices of the triangles in the neighbours.
  for( i = 0; i < surface->n_items; i++ ) {
    for( j = 0; j < 3; j++ ) {
      jj = surface->indices[3*i+j];
      for( k = 0; k < n_ngh[jj]; k++ ) {
        if( ngh[jj][k] == -1 ) {
          ngh[jj][k] = i;
          break;
        }
      }
    }
  }

  // Now create a sort closed loop of the node neighbours.
  // This is needed by the parametric=0 FEM algorithm.
  //
  //         1 ----- 2
  //          /\   /\
  //         /  \ /  \
  //       0 ----P---- 3
  //         \  / \  /
  //          \/   \/
  //         5 ----- 4
  //

  int * tmp;
  ALLOC( tmp, 2*max_ngh );

  for( i = 0; i < surface->n_points; i++ ) {
    for( k = 0; k < n_ngh[i]; k++ ) {
      tri = &(surface->indices[3*ngh[i][k]]);
      for( j = 0; j < 3; j++ ) {
        if( tri[j] == i ) break;
      }
      tmp[2*k+0] = tri[(j+1)%3];
      tmp[2*k+1] = tri[(j+2)%3];
    }

    ngh[i][0] = tmp[0];
    ngh[i][1] = tmp[1];
    for( k = 2; k < n_ngh[i]; k++ ) {
      for( j = 1; j < n_ngh[i]; j++ ) {
        if( tmp[2*j] == ngh[i][k-1] || tmp[2*j+1] == ngh[i][k-1] ) {
          if( tmp[2*j] == ngh[i][k-1] ) {
            ngh[i][k] = tmp[2*j+1];
          } else {
            ngh[i][k] = tmp[2*j];
          }
          tmp[2*j] = -1;
          tmp[2*j+1] = -1;
          break;
        }
      }
    }
  }

  *n_neighbours_return = n_ngh;
  *neighbours_return = ngh;

  FREE( tmp );

  return OK;

}

// -------------------------------------------------------------------
// Read a scalar from a file.
//
static Real * read_scalar( int n_points, char * in_file ) {

  int i;

  Real * scalar = (Real*)malloc( n_points * sizeof( Real ) );
  FILE * fp = fopen( in_file, "r" );
  for( i = 0; i < n_points; i++ ) {
    fscanf( fp, "%lg\n", &scalar[i] );
  }
  fclose( fp );
  return( scalar );
}

// -------------------------------------------------------------------
// Save a scalar to a file.
//
private void save_scalar( int n_points, Real scalar[], char * out_file ) {

  int i;

  FILE * fp = fopen( out_file, "w" );
  for( i = 0; i < n_points; i++ ) {
    fprintf( fp, "%g\n", scalar[i] );
  }
  fclose( fp );
}

