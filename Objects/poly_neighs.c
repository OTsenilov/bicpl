
#include  <def_mni.h>
#include  <def_alloc.h>
#include  <def_objects.h>
#include  <def_progress.h>
#include  <def_hash.h>

#define  INVALID_ID       -1

typedef  struct
{
    int   polygon_index;   
} edge_info_struct;

#define  INITIAL_HASH_TABLE_SIZE   2.0   /* times number of polygons */
#define  ENLARGE_THRESHOLD         0.25
#define  NEW_DENSITY               0.125

static    void                assign_neighbours();

public  Status  check_polygons_neighbours_computed( polygons )
    polygons_struct   *polygons;
{
    Status  status;
    Status  create_polygon_neighbours();

    status = OK;

    if( polygons->neighbours == (int *) 0 )
        status = create_polygon_neighbours( polygons->n_items,
                                            polygons->indices,
                                            polygons->end_indices,
                                            &polygons->neighbours );

    return( status );
}

private   Status   create_polygon_neighbours( n_polygons, indices, end_indices,
                                              neighbours )
    int    n_polygons;
    int    indices[];
    int    end_indices[];
    int    *neighbours[];
{
    int                 keys[2], i, edge, i0, i1, size;
    int                 start_index, end_index;
    edge_info_struct    *edge_ptr;
    hash_table_struct   edge_table;
    hash_table_pointer  hash_ptr;
    Status              status;
    progress_struct     progress;

    if( n_polygons > 0 )
    {
        ALLOC( status, *neighbours, end_indices[n_polygons-1] );

        if( status == OK )
        {
            for_less( i, 0, end_indices[n_polygons-1] )
            {
                (*neighbours)[i] = INVALID_ID;
            }
        }
    }

    if( status == OK )
    {
        status = initialize_hash_table( &edge_table, 2,
                                    (int) (INITIAL_HASH_TABLE_SIZE*n_polygons),
                                    ENLARGE_THRESHOLD, NEW_DENSITY );
    }

    if( status == OK )
    {
        end_index = 0;

        initialize_progress_report( &progress, FALSE, n_polygons,
                                    "Neighbour-finding" );

        for_less( i, 0, n_polygons )
        {
            start_index = end_index;
            end_index = end_indices[i];

            size = end_index - start_index;

            for_less( edge, 0, size )
            {
                i0 = indices[start_index+edge];
                i1 = indices[start_index+(edge+1) % size];

                keys[0] = MIN( i0, i1 );
                keys[1] = MAX( i0, i1 );

                if( remove_from_hash_table( &edge_table, keys,
                                            (void **) &edge_ptr ) )
                {
                    assign_neighbours( indices, end_indices, *neighbours,
                                       i, edge, edge_ptr->polygon_index );
                    FREE( status, edge_ptr );
                }
                else
                {
                    ALLOC( status, edge_ptr, 1 );
                    if( status == OK )
                    {
                        edge_ptr->polygon_index = i;
                        status = insert_in_hash_table( &edge_table, keys,
                                                       (void *) edge_ptr );
                    }
                }

                if( status != OK )
                {
                    break;
                }
            }

            if( status != OK )
            {
                break;
            }

            update_progress_report( &progress, i+1 );
        }

        terminate_progress_report( &progress );
    }

    if( status == OK )
    {
        end_index = 0;

        initialize_hash_pointer( &hash_ptr );

        while( status == OK && get_next_hash_entry( &edge_table, &hash_ptr,
                                            (void **) &edge_ptr ) )
        {
            FREE( status, edge_ptr );
        }
    }

    if( status == OK )
    {
        status = delete_hash_table( &edge_table );
    }

    return( status );
}

public   Status   free_polygon_neighbours( neighbours )
    int    neighbours[];
{
    Status   status;

    FREE( status, neighbours );

    return( status );
}

private  void  assign_neighbours( indices, end_indices, neighbours,
                                  polygon1, edge1, polygon2 )
    int       indices[];
    int       end_indices[];
    int       neighbours[];
    int       polygon1;
    int       edge1;
    int       polygon2;
{
    int   i0, i1, edge2, start_index1, size1, start_index2, size2;
    int   p2_i0, p2_i1;

    if( polygon1 == 0 )
        start_index1 = 0;
    else
        start_index1 = end_indices[polygon1-1];
    size1 = end_indices[polygon1] - start_index1;

    if( polygon2 == 0 )
        start_index2 = 0;
    else
        start_index2 = end_indices[polygon2-1];
    size2 = end_indices[polygon2] - start_index2;

    i0 = indices[start_index1 + edge1];
    i1 = indices[start_index1 + (edge1 + 1) % size1];

    for_less( edge2, 0, size2 )
    {
        p2_i0 = indices[start_index2 + edge2];
        p2_i1 = indices[start_index2 + (edge2 + 1) % size2];

        if( (p2_i0 == i0 && p2_i1 == i1) ||
            (p2_i0 == i1 && p2_i1 == i0) )
        {
            break;
        }
    }

    if( edge2 == size2 )
    {
        HANDLE_INTERNAL_ERROR( "assign neighbours" );
    }

    neighbours[POINT_INDEX( end_indices, polygon1, edge1 )] = polygon2;
    neighbours[POINT_INDEX( end_indices, polygon2, edge2 )] = polygon1;
}