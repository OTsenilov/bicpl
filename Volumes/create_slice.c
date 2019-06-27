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

#include  "bicpl_internal.h"

#ifndef lint
static char rcsid[] = "$Header: /home/users/clepage/CVS/libraries/bicpl/Volumes/create_slice.c,v 1.47 2005/08/17 22:26:19 bert Exp $";
#endif

static void  create_pixel_mapping(
    Volume          volume1,
    int             n_slices1,
    Real            **origins1,
    Real            x_axis1[],
    Real            y_axis1[],
    Real            x_translation1,
    Real            y_translation1,
    Real            x_scale1,
    Real            y_scale1,
    Volume          volume2,
    int             n_slices2,
    Real            **origins2,
    Real            x_axis2[],
    Real            y_axis2[],
    Real            x_translation2,
    Real            y_translation2,
    Real            x_scale2,
    Real            y_scale2,
    Real            real_x_axis1[],
    Real            real_y_axis1[],
    Real            ***real_origins1,
    Real            real_x_axis2[],
    Real            real_y_axis2[],
    Real            ***real_origins2 )
{
    int          n_dimensions1, n_dimensions2, s;

    n_dimensions1 = get_volume_n_dimensions( volume1 );

    ALLOC2D( *real_origins1, n_slices1, n_dimensions1 );

    for_less( s, 0, n_slices1 )
    {
        get_mapping( volume1, origins1[s], x_axis1, y_axis1,
                     x_translation1, y_translation1, x_scale1, y_scale1,
                     (*real_origins1)[s], real_x_axis1, real_y_axis1 );
    }

    if( volume2 != NULL )
    {
        n_dimensions2 = get_volume_n_dimensions( volume2 );

        ALLOC2D( *real_origins2, n_slices2, n_dimensions2 );

        for_less( s, 0, n_slices2 )
        {
            get_mapping( volume2, origins2[s], x_axis2, y_axis2,
                         x_translation2, y_translation2, x_scale2, y_scale2,
                         (*real_origins2)[s], real_x_axis2, real_y_axis2 );
        }
    }
}

static void  set_pixel_range(
    Volume          volume1,
    int             n_slices1,
    Real            **real_origins1,
    Real            real_x_axis1[],
    Real            real_y_axis1[],
    Volume          volume2,
    int             n_slices2,
    Real            **real_origins2,
    Real            real_x_axis2[],
    Real            real_y_axis2[],
    int             x_viewport_size,
    int             y_viewport_size,
    Pixel_types     pixel_type,
    int             *n_pixels_alloced,
    pixels_struct   *pixels )
{
    int          s, x_size, y_size;
    int          x_pixel_start, x_pixel_end, y_pixel_start, y_pixel_end;

    x_pixel_start = 0;
    x_pixel_end = x_viewport_size - 1;
    y_pixel_start = 0;
    y_pixel_end = y_viewport_size - 1;

    for_less( s, 0, n_slices1 )
    {
        clip_viewport_to_volume( volume1, real_origins1[s],
                                 real_x_axis1, real_y_axis1,
                                 &x_pixel_start, &x_pixel_end,
                                 &y_pixel_start, &y_pixel_end );
    }

    if( volume2 != NULL )
    {
        for_less( s, 0, n_slices2 )
        {
            clip_viewport_to_volume( volume2, real_origins2[s],
                                     real_x_axis2, real_y_axis2,
                                     &x_pixel_start, &x_pixel_end,
                                     &y_pixel_start, &y_pixel_end );
        }
    }

    x_size = x_pixel_end - x_pixel_start + 1;
    y_size = y_pixel_end - y_pixel_start + 1;
    if( x_size < 0 )
        x_size = 0;
    if( y_size < 0 )
        y_size = 0;

    modify_pixels_size( n_pixels_alloced, pixels, x_size, y_size, pixel_type );
    pixels->x_position = x_pixel_start;
    pixels->y_position = y_pixel_start;
    pixels->x_zoom = 1.0;
    pixels->y_zoom = 1.0;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_weighted_volume_slices
@INPUT      : volume1,
              n_slices1   - number of parallel slices to weighted-sum.
              origins1
              x_axis1
              y_axis1
              weights1
              x_translation1
              y_translation1
              x_scale1
              y_scale1
              volume2,
              n_slices2   - number of parallel slices to weighted-sum.
              origins2
              x_axis2
              y_axis2
              weights2
              x_translation2
              y_translation2
              x_scale2
              y_scale2
              x_viewport_size
              y_viewport_size
              pixel_type
              degrees_continuity
              cmode_colour_map
              rgb_colour_map
              empty_colour
              n_pixels_alloced
@OUTPUT     : pixels
@RETURNS    : 
@DESCRIPTION: Renders a volume or a merge of 2 volumes to pixels.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static void  create_weighted_volume_slices(
    Volume          volume1,
    int             n_slices1,
    Real            **origins1,
    Real            x_axis1[],
    Real            y_axis1[],
    Real            weights1[],
    Volume          volume2,
    int             n_slices2,
    Real            **origins2,
    Real            x_axis2[],
    Real            y_axis2[],
    Real            weights2[],
    int             x_pixel_start,
    int             x_pixel_end,
    int             y_pixel_start,
    int             y_pixel_end,
    int             degrees_continuity,
    unsigned short  **cmode_colour_map,
    Colour          **rgb_colour_map,
    Colour          empty_colour,
    void            *render_storage,
    pixels_struct   *pixels )
{
    Data_types   volume1_data_type, volume2_data_type;
    void         *volume_data1, *volume_data2;
    int          strides1[MAX_DIMENSIONS], strides2[MAX_DIMENSIONS];
    int          n_dimensions1, n_dimensions2;
    int          axis, s, c;
    int          sizes1[MAX_DIMENSIONS], sizes2[MAX_DIMENSIONS];

    if( pixels->x_size > 0 && pixels->y_size > 0 )
    {
        n_dimensions1 = get_volume_n_dimensions( volume1 );

        for_less( s, 0, n_slices1 )
        {
            for_less( c, 0, n_dimensions1 )
            {
                origins1[s][c] += (Real) pixels->x_position * x_axis1[c] +
                                  (Real) pixels->y_position * y_axis1[c];
            }
        }

        if( !volume1->is_cached_volume ) {
            GET_VOXEL_PTR( volume_data1, volume1, 0, 0, 0, 0, 0 );
	}
        volume1_data_type = get_volume_data_type( volume1 );

        get_volume_sizes( volume1, sizes1 );

        strides1[n_dimensions1-1] = 1;
        for_down( axis, n_dimensions1 - 2,  0 )
             strides1[axis] = strides1[axis+1] * sizes1[axis+1];

        if( volume2 != NULL )
        {
            n_dimensions2 = get_volume_n_dimensions( volume2 );

            for_less( s, 0, n_slices2 )
            {
                for_less( c, 0, n_dimensions2 )
                {
                    origins2[s][c] += (Real) pixels->x_position * x_axis2[c] +
                                      (Real) pixels->y_position * y_axis2[c];
                }
            }

            if( !volume2->is_cached_volume ) {
                GET_VOXEL_PTR( volume_data2, volume2, 0, 0, 0, 0, 0 );
	    }
            volume2_data_type = get_volume_data_type( volume2 );
            get_volume_sizes( volume2, sizes2 );

            strides2[n_dimensions2-1] = 1;
            for_down( axis, n_dimensions2 - 2, 0 )
                 strides2[axis] = strides2[axis+1] * sizes2[axis+1];
        }
        else
            volume_data2 = NULL;

        if( x_pixel_start > x_pixel_end || y_pixel_start > y_pixel_end )
        {
            x_pixel_start = 0;
            x_pixel_end = pixels->x_size-1;
            y_pixel_start = 0;
            y_pixel_end = pixels->y_size-1;
        }
        else
        {
            if( x_pixel_start < 0 )
                x_pixel_start = 0;
            if( x_pixel_end >= pixels->x_size )
                x_pixel_end = pixels->x_size-1;

            if( y_pixel_start < 0 )
                y_pixel_start = 0;
            if( y_pixel_end >= pixels->y_size )
                y_pixel_end = pixels->y_size-1;
        }

        if( volume1->is_cached_volume ||
            (volume2 != NULL && volume2->is_cached_volume) ||
            (degrees_continuity != -1 && n_slices1 == 1 &&
            (volume_data2 == NULL || n_slices2 == 1)) )
        {
            interpolate_volume_to_slice( volume1, n_dimensions1,
                                         origins1[0],
                                         x_axis1, y_axis1,
                                         volume2, n_dimensions2,
                                         (volume2 == NULL) ? NULL : origins2[0],
                                         x_axis2, y_axis2,
                                         x_pixel_start, x_pixel_end,
                                         y_pixel_start, y_pixel_end,
                                         degrees_continuity,
                                         cmode_colour_map, rgb_colour_map,
                                         empty_colour, pixels );
        }
        else
        {
            render_volume_to_slice( n_dimensions1, sizes1, volume_data1,
                                    volume1_data_type, n_slices1, weights1,
                                    strides1, origins1, x_axis1, y_axis1,
                                    n_dimensions2, sizes2, volume_data2,
                                    volume2_data_type, n_slices2, weights2,
                                    strides2, origins2, x_axis2, y_axis2,
                                    x_pixel_start, x_pixel_end,
                                    y_pixel_start, y_pixel_end,
                                    cmode_colour_map,
                                    rgb_colour_map, empty_colour,
                                    render_storage, pixels );
        }
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_filter_slices
@INPUT      : volume
              position
              x_axis
              y_axis
              filter_type
              filter_width
@OUTPUT     : n_slices
              origins
              weights
@RETURNS    : TRUE if successful
@DESCRIPTION: Converts the filter type and width to a list of parallel slices
              to be used in a weighted-sum rendering.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static BOOLEAN  get_filter_slices(
    Volume          volume,
    Real            position[],
    Real            x_axis[],
    Real            y_axis[],
    Filter_types    filter_type,
    Real            filter_width,
    int             *n_slices,
    Real            ***origins,
    Real            **weights )
{
    int    c, a1, a2;
    Real   direction[N_DIMENSIONS], separations[MAX_DIMENSIONS];

    if( filter_type != NEAREST_NEIGHBOUR )
    {
        if( get_volume_n_dimensions( volume ) != N_DIMENSIONS )
        {
            print_error(
              "If volume not 3D, can only do nearest neighbour filtering.\n" );
            return( FALSE );
        }

        get_volume_separations( volume, separations );
        for_less( c, 0, N_DIMENSIONS )
        {
            a1 = (c + 1) % N_DIMENSIONS;
            a2 = (c + 2) % N_DIMENSIONS;
            direction[c] = x_axis[a1] * y_axis[a2] - x_axis[a2] * y_axis[a1];
            direction[c] *= FABS( separations[a1] * separations[a2] /
                                  separations[c] );
        }
    }

    *n_slices = get_slice_weights_for_filter( volume, position, direction,
                                              filter_type, filter_width,
                                              origins, weights );

    return( TRUE );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_volume_slice
@INPUT      : volume1          - the volume to create a slice for
              filter_type1     - filter_type, usually NEAREST_NEIGHBOUR
              filter_width1    - width of filter, used for BOX, TRIANGLE, GAUSS
              slice_position1  - the voxel coordinate of the slice
              x_axis1          - the x axis in voxels
              y_axis1          - the y axis in voxels
              x_translation1   - pixel translation for viewing
              y_translation1   - pixel translation for viewing
              x_scale1         - pixel zoom for viewing
              y_scale1         - pixel zoom for viewing
              volume2          - second volume to be merged with first, or null
              filter_type2     - filter_type, usually NEAREST_NEIGHBOUR
              filter_width2    - width of filter, used for BOX, TRIANGLE, GAUSS
              slice_position2  - the voxel coordinate of the slice
              x_axis2          - the x axis in voxels
              y_axis2          - the y axis in voxels
              x_translation2   - pixel translation for viewing
              y_translation2   - pixel translation for viewing
              x_scale2         - pixel zoom for viewing
              y_scale2         - pixel zoom for viewing
              x_axis_index     - X,Y, or Z
              y_axis_index     - X,Y, or Z
              axis_index       - X,Y, or Z
              x_viewport_size  - will be clipped to this size
              y_viewport_size  - will be clipped to this size
              pixel_type       - RGB_PIXEL or COLOUR_INDEX_PIXEL for rgb/cmap
              interpolation_flag - ignored for now
              cmode_colour_map - if pixel_type == COLOUR_INDEX_PIXEL, then
                          2d array of 16 bit colour indices for merged slices,
                          or pointer to 1d array of colour indices for volume1
              rgb_colour_map - if pixel_type == RGB_PIXEL, then
                          2d array of 24 bit colours for merged slices,
                          or pointer to 1d array of colours for volume1
@OUTPUT     : n_pixels_alloced - a pointer to the size alloced.  Before first
                          call, set size alloced to zero, and all calls,
                          pass pointer to size alloced, and pointer to pixels.
              pixels           - 2d pixels array created, and realloced as
                                 necessary, assuming, n_pixels_alloced is a
                                 pointer to the current alloc size of pixels.
@RETURNS    : 
@DESCRIPTION: Creates a slice of one volume or merged slice of two, suitable
              for graphics display.
@CREATED    : Mar   1993           David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI void  create_volume_slice(
    Volume          volume1,
    Filter_types    filter_type1,
    Real            filter_width1,
    Real            slice_position1[],
    Real            x_axis1[],
    Real            y_axis1[],
    Real            x_translation1,
    Real            y_translation1,
    Real            x_scale1,
    Real            y_scale1,
    Volume          volume2,
    Filter_types    filter_type2,
    Real            filter_width2,
    Real            slice_position2[],
    Real            x_axis2[],
    Real            y_axis2[],
    Real            x_translation2,
    Real            y_translation2,
    Real            x_scale2,
    Real            y_scale2,
    int             x_viewport_size,
    int             y_viewport_size,
    int             x_pixel_start,
    int             x_pixel_end,
    int             y_pixel_start,
    int             y_pixel_end,
    Pixel_types     pixel_type,
    int             degrees_continuity,
    unsigned short  **cmode_colour_map,
    Colour          **rgb_colour_map,
    Colour          empty_colour,
    void            *render_storage,
    BOOLEAN         clip_pixels_flag,
    int             *n_pixels_alloced,
    pixels_struct   *pixels )
{
    int          n_slices1, n_slices2;
    Real         **positions1, **positions2, *weights1, *weights2;
    Real         **real_origins1, **real_origins2;
    Real         real_x_axis1[MAX_DIMENSIONS], real_y_axis1[MAX_DIMENSIONS];
    Real         real_x_axis2[MAX_DIMENSIONS], real_y_axis2[MAX_DIMENSIONS];

    real_origins1 = NULL;       /* Initialize this!! */
    real_origins2 = NULL;       /* Initialize this!! */

    if( !get_filter_slices( volume1, slice_position1, x_axis1, y_axis1,
                            filter_type1, filter_width1, &n_slices1,
                            &positions1, &weights1 ) )
    {
        modify_pixels_size( n_pixels_alloced, pixels, 0, 0, pixel_type );
        return;
    }

    if( volume2 != NULL )
    {
        if( !get_filter_slices( volume2, slice_position2, x_axis2, y_axis2,
                                filter_type2, filter_width2, &n_slices2,
                                &positions2, &weights2 ) )
        {
            modify_pixels_size( n_pixels_alloced, pixels, 0, 0, pixel_type );
            return;
        }
    }

    create_pixel_mapping( volume1, n_slices1, positions1, x_axis1, y_axis1,
                          x_translation1, y_translation1, x_scale1, y_scale1,
                          volume2, n_slices2, positions2, x_axis2, y_axis2,
                          x_translation2, y_translation2, x_scale2, y_scale2,
                          real_x_axis1, real_y_axis1, &real_origins1,
                          real_x_axis2, real_y_axis2, &real_origins2 );

    if( clip_pixels_flag )
    {
        set_pixel_range( volume1, n_slices1,
                         real_origins1, real_x_axis1, real_y_axis1,
                         volume2, n_slices2,
                         real_origins2, real_x_axis2, real_y_axis2,
                         x_viewport_size, y_viewport_size,
                         pixel_type, n_pixels_alloced, pixels );
    }

    create_weighted_volume_slices( volume1, n_slices1,
                                   real_origins1, real_x_axis1, real_y_axis1,
                                   weights1,
                                   volume2, n_slices2,
                                   real_origins2, real_x_axis2, real_y_axis2,
                                   weights2,
                                   x_pixel_start, x_pixel_end,
                                   y_pixel_start, y_pixel_end,
                                   degrees_continuity,
                                   cmode_colour_map, rgb_colour_map,
                                   empty_colour, render_storage,
                                   pixels );

    if( volume2 != NULL )
    {
        FREE2D( positions2 );
        FREE( weights2 );
        FREE2D( real_origins2 );
    }

    FREE2D( positions1 );
    FREE( weights1 );
    FREE2D( real_origins1 );
}

BICAPI void  set_volume_slice_pixel_range(
    Volume          volume1,
    Filter_types    filter_type1,
    Real            filter_width1,
    Real            slice_position1[],
    Real            x_axis1[],
    Real            y_axis1[],
    Real            x_translation1,
    Real            y_translation1,
    Real            x_scale1,
    Real            y_scale1,
    Volume          volume2,
    Filter_types    filter_type2,
    Real            filter_width2,
    Real            slice_position2[],
    Real            x_axis2[],
    Real            y_axis2[],
    Real            x_translation2,
    Real            y_translation2,
    Real            x_scale2,
    Real            y_scale2,
    int             x_viewport_size,
    int             y_viewport_size,
    Pixel_types     pixel_type,
    int             *n_pixels_alloced,
    pixels_struct   *pixels )
{
    int          n_slices1, n_slices2;
    Real         **positions1, **positions2, *weights1, *weights2;
    Real         **real_origins1, **real_origins2;
    Real         real_x_axis1[MAX_DIMENSIONS], real_y_axis1[MAX_DIMENSIONS];
    Real         real_x_axis2[MAX_DIMENSIONS], real_y_axis2[MAX_DIMENSIONS];

    if( !get_filter_slices( volume1, slice_position1, x_axis1, y_axis1,
                            filter_type1, filter_width1, &n_slices1,
                            &positions1, &weights1 ) )
    {
        modify_pixels_size( n_pixels_alloced, pixels, 0, 0, pixel_type );
        return;
    }

    if( volume2 != NULL )
    {
        if( !get_filter_slices( volume2, slice_position2, x_axis2, y_axis2,
                                filter_type2, filter_width2, &n_slices2,
                                &positions2, &weights2 ) )
        {
            modify_pixels_size( n_pixels_alloced, pixels, 0, 0, pixel_type );
            return;
        }
    }

    create_pixel_mapping( volume1, n_slices1, positions1, x_axis1, y_axis1,
                          x_translation1, y_translation1, x_scale1, y_scale1,
                          volume2, n_slices2, positions2, x_axis2, y_axis2,
                          x_translation2, y_translation2, x_scale2, y_scale2,
                          real_x_axis1, real_y_axis1, &real_origins1,
                          real_x_axis2, real_y_axis2, &real_origins2 );

    set_pixel_range( volume1, n_slices1,
                     real_origins1, real_x_axis1, real_y_axis1,
                     volume2, n_slices2,
                     real_origins2, real_x_axis2, real_y_axis2,
                     x_viewport_size, y_viewport_size,
                     pixel_type, n_pixels_alloced, pixels );

    if( volume2 != NULL )
    {
        FREE2D( positions2 );
        FREE( weights2 );
        FREE2D( real_origins2 );
    }

    FREE2D( positions1 );
    FREE( weights1 );
    FREE2D( real_origins1 );
}
