#include  <module.h>

public  void   create_slice_quadmesh(
    Volume           volume,
    int              axis_index,
    Real             voxel_position,
    quadmesh_struct  *quadmesh )
{
    Point            point;
    Vector           normal;
    Real             x_w, y_w, z_w;
    int              sizes[N_DIMENSIONS];
    Real             voxel[N_DIMENSIONS];
    int              x, y, x_axis, y_axis;
    Surfprop         spr;

    x_axis = (axis_index + 1) % N_DIMENSIONS;
    y_axis = (axis_index + 2) % N_DIMENSIONS;
    get_volume_sizes( volume, sizes );

    get_default_surfprop( &spr );
    initialize_quadmesh( quadmesh, WHITE, &spr, sizes[x_axis], sizes[y_axis] );

    voxel[axis_index] = voxel_position;

    fill_Vector( normal, 0.0, 0.0, 0.0 );
    Vector_coord( normal, axis_index ) = 1.0;

    for_less( x, 0, sizes[x_axis] )
    {
        voxel[x_axis] = (Real) x;
        for_less( y, 0, sizes[y_axis] )
        {
            voxel[y_axis] = (Real) y;

            convert_voxel_to_world( volume, voxel, &x_w, &y_w, &z_w );
            fill_Point( point, x_w, y_w, z_w );
            set_quadmesh_point( quadmesh, x, y, &point, &normal );
        }
    }
}
