
#include  <def_mni.h>

#define  TALAIRACH_OFFSET   16.0

#define  FALSE  0
#define  TRUE   1

private  int   nx = 128;
private  int   ny = 128;
private  int   nz = 80;

private  Real  x_low  = -1.28;
private  Real  x_high =  1.28;
private  Real  y_low  = -1.28;
private  Real  y_high =  1.28;
private  Real  z_low  = -0.50;
private  Real  z_high =  1.10;

private  Real  x_dist_minus_1 = -67.0;
private  Real  x_dist_1       =  67.0;
private  Real  y_dist_minus_1 = -86.0;
private  Real  y_dist_1       =  86.0;
private  Real  z_dist_0       =   0.0;
private  Real  z_dist_1       =  75.0;

private  Boolean   initialized = FALSE;

private  void  read_talairach_coordinate_system( void );

private  void  check_initialized( void )
{
    if( !initialized )
    {
        read_talairach_coordinate_system();
        initialized = TRUE;
    }
}

private  void  read_talairach_coordinate_system( void )
{
    int     i, n_lines, ch;
    FILE    *file;
    char    *name = "brain_size.dat";
    char    filename[234];
    int     okay;
    int     tmp_nx, tmp_ny, tmp_nz;
    Real    tmp_x_low, tmp_x_high, tmp_y_low, tmp_y_high, tmp_z_low;
    Real    tmp_z_high, tmp_x_dist_minus_1, tmp_x_dist_1, tmp_y_dist_minus_1;
    Real    tmp_y_dist_1, tmp_z_dist_0, tmp_z_dist_1;

    okay = TRUE;

    if( file_exists(name) )
    {
        (void) strcpy( filename, name );
    }
    else
    {
        if( getenv( "PET_ROI" ) == (char *) 0 )
        {
            print( "You have to setenv PET_ROI.\n" );
            okay = FALSE;
        }
        else
        {
            get_absolute_filename( name, getenv("PET_ROI"), filename );
        }
    }

    if( okay )
    {
        file = fopen( filename, "r" );

        if( file == (FILE *) 0 )
        {
            print( "Cannot open %s\n", filename );
            okay = FALSE;
        }
    }

    if( okay )
    {
        (void) input_int( file, &n_lines );

        for( i = 0;  i < n_lines+1;  ++i )
        {
            while( (ch = fgetc(file)) != EOF && ch != '\n' )
            {}
        }

        if( input_int( file, &tmp_nx ) != OK ||
            input_int( file, &tmp_ny ) != OK ||
            input_int( file, &tmp_nz ) != OK ||
            input_real( file, &tmp_x_low ) != OK ||
            input_real( file, &tmp_x_high ) != OK ||
            input_real( file, &tmp_y_low ) != OK ||
            input_real( file, &tmp_y_high ) != OK ||
            input_real( file, &tmp_z_low ) != OK ||
            input_real( file, &tmp_z_high ) != OK ||
            input_real( file, &tmp_x_dist_minus_1 ) != OK ||
            input_real( file, &tmp_x_dist_1 ) != OK ||
            input_real( file, &tmp_y_dist_minus_1 ) != OK ||
            input_real( file, &tmp_y_dist_1 ) != OK ||
            input_real( file, &tmp_z_dist_0 ) != OK ||
            input_real( file, &tmp_z_dist_1 ) != OK )
        {
            print( "Error reading info from %s\n", filename );
            okay = FALSE;
        }
    }

    if( okay )
    {
        nx = tmp_nx;
        ny = tmp_ny;
        nz = tmp_nz;

        x_low = tmp_x_low;
        x_high = tmp_x_high;
        y_low = tmp_y_low;
        y_high = tmp_y_high;
        z_low = tmp_z_low;
        z_high = tmp_z_high;

        x_dist_minus_1 = tmp_x_dist_minus_1;
        x_dist_1 = tmp_x_dist_1;
        y_dist_minus_1 = tmp_y_dist_minus_1;
        y_dist_1 = tmp_y_dist_1;
        z_dist_0 = tmp_z_dist_0;
        z_dist_1 = tmp_z_dist_1;

        (void) fclose( file );
    }
}

public  void  convert_voxel_to_talairach(
    Real   x_voxel,
    Real   y_voxel,
    Real   z_voxel,
    int    nx_voxels,
    int    ny_voxels,
    int    nz_voxels,
    Real   *x_tal,
    Real   *y_tal,
    Real   *z_tal )
{
    check_initialized();

    *x_tal = x_voxel / (Real) nx_voxels * (Real) nx;
    *y_tal = y_voxel / (Real) ny_voxels * (Real) ny;
    *z_tal = z_voxel / (Real) nz_voxels * (Real) nz;
}

public  void  convert_talairach_to_voxel(
    Real   x_tal,
    Real   y_tal,
    Real   z_tal,
    int    nx_voxels,
    int    ny_voxels,
    int    nz_voxels,
    Real   *x_voxel,
    Real   *y_voxel,
    Real   *z_voxel )
{
    check_initialized();

    *x_voxel = x_tal / (Real) nx * (Real) nx_voxels;
    *y_voxel = y_tal / (Real) ny * (Real) ny_voxels;
    *z_voxel = z_tal / (Real) nz * (Real) nz_voxels;
}

private  Real  convert_to_mm(
    Real    tal,
    int     n_planes,
    Real    limit_low,
    Real    limit_high,
    Real    brain_dist_low,
    Real    brain_dist_high,
    Real    brain_limit_low,
    Real    brain_limit_high,
    Real    *mm )
{
    Real   stereotactic;

    stereotactic = limit_low + (limit_high - limit_low) * tal / (Real) n_planes;

    *mm = brain_dist_low + (brain_dist_high - brain_dist_low) *
      (stereotactic - brain_limit_low) / (brain_limit_high - brain_limit_low);
}

public  void  convert_talairach_to_mm(
    Real   x_tal,
    Real   y_tal,
    Real   z_tal,
    Real   *x_mm,
    Real   *y_mm,
    Real   *z_mm )
{
    check_initialized();

    convert_to_mm( x_tal, nx, x_low, x_high, x_dist_minus_1, x_dist_1,
                   -1.0, 1.0, x_mm );
    convert_to_mm( y_tal, ny, y_low, y_high, y_dist_minus_1, y_dist_1,
                   -1.0, 1.0, y_mm );
    *y_mm -= TALAIRACH_OFFSET;
    convert_to_mm( z_tal, nz, z_low, z_high, z_dist_0, z_dist_1,
                   0.0, 1.0, z_mm );
}

private  Real  convert_from_mm(
    Real    mm,
    int     n_planes,
    Real    limit_low,
    Real    limit_high,
    Real    brain_dist_low,
    Real    brain_dist_high,
    Real    brain_limit_low,
    Real    brain_limit_high,
    Real    *tal )
{
    Real  stereotactic;

    stereotactic = brain_limit_low +
                   (mm - brain_dist_low) * (brain_limit_high - brain_limit_low)
                   / (brain_dist_high - brain_dist_low);

    *tal = (stereotactic - limit_low) * (Real) n_planes /
           (limit_high - limit_low);
}

public  void  convert_mm_to_talairach(
    Real   x_mm,
    Real   y_mm,
    Real   z_mm,
    Real   *x_tal,
    Real   *y_tal,
    Real   *z_tal )
{
    check_initialized();

    convert_from_mm( x_mm, nx, x_low, x_high, x_dist_minus_1, x_dist_1,
                     -1.0, 1.0, x_tal );
    convert_from_mm( y_mm + TALAIRACH_OFFSET, ny, y_low, y_high,
                     y_dist_minus_1, y_dist_1, -1.0, 1.0, y_tal );
    convert_from_mm( z_mm, nz, z_low, z_high, z_dist_0, z_dist_1,
                     0.0, 1.0, z_tal );
}