
/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      :  
@CREATED    : Sep. 10, 1996    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */


#ifndef lint
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Geometry/solve_plane.c,v 1.10 1996-09-18 18:14:42 david Exp $";
#endif

#include  <internal_volume_io.h>
#include  <geom.h>
#include  <trans.h>

public  BOOLEAN  get_interpolation_weights_2d(
    Real   x,
    Real   y,
    int    n_points,
    Real   xs[],
    Real   ys[],
    Real   weights[] )
{
    int   i;
    Real  x_fact, y_fact, constant, n;
    Real  aa, ab, ac, bb, bc, dx, dy;
    Real  acbb, abac, aabb, denom, aabc, abab;

    aa = 0.0;
    ab = 0.0;
    ac = 0.0;
    bb = 0.0;
    bc = 0.0;

    for_less( i, 0, n_points )
    {
        dx = xs[i] - x;
        dy = ys[i] - y;
        aa += dx * dx;
        ab += dx * dy;
        ac += dx;
        bb += dy * dy;
        bc += dy;
    }

    n = (Real) n_points;
    aabb = aa * bb;
    acbb = ac * bb;
    aabc = aa * bc;
    abab = ab * ab;
    abac = ab * ac;

    denom = -aabb * n + ac * acbb + bc * aabc + abab*n - 2.0 * abac * bc;

    if( denom == 0.0 )
        return( FALSE );

    x_fact = (acbb - ab * bc) / denom;
    y_fact = (aabc - abac) / denom;
    constant = (abab - aabb) / denom - x * x_fact - y * y_fact;

    for_less( i, 0, n_points )
        weights[i] = constant + y_fact * ys[i] + x_fact * xs[i];

    return( TRUE );
}

#ifdef   DEBUG
#define  DEBUG

#include <prog_utils.h>
#include <numerical.h>

private  void  test_solution(
    int    dim,
    Real   x,
    Real   y,
    int    n_points,
    Real   xs[],
    Real   ys[],
    Real   x_weights[],
    Real   y_weights[],
    Real   constant )
{
    int        iter, n_iters, p;
    Real       angle, x_trans, y_trans, xt, yt, zt, correct, value;
    Transform  transform, rotation, translation;

    n_iters = 100;

    for_less( iter, 0, n_iters )
    {
        angle = 2.0 * PI * get_random_0_to_1();
        x_trans = 10.0 * get_random_0_to_1() - 5.0;
        y_trans = 10.0 * get_random_0_to_1() - 5.0;

        make_rotation_transform( angle, Z, &rotation );
        make_translation_transform( x_trans, y_trans, 0.0, &translation );
        concat_transforms( &transform, &translation, &rotation );

        value = constant;

        for_less( p, 0, n_points )
        {
            transform_point( &transform, xs[p], ys[p], 0.0, &xt, &yt, &zt );
            value += x_weights[p] * xt + y_weights[p] * yt;
        }

        transform_point( &transform, x, y, 0.0, &xt, &yt, &zt );
        correct = (dim == 0) ? xt : yt;

        if( !numerically_close( value, correct, 1.0e-6 ) )
        {
            print( "get_prediction_weights_2d_for_1_coord: %g %g\n",
                   correct, value );
            break;
        }
    }
}
#endif

private  BOOLEAN   get_two_point_prediction(
    Real   x,
    Real   y,
    Real   x1,
    Real   y1,
    Real   x2,
    Real   y2,
    Real   *xwx1,
    Real   *xwy1,
    Real   *xwx2,
    Real   *xwy2,
    Real   *ywx1,
    Real   *ywy1,
    Real   *ywx2,
    Real   *ywy2 )
{
    Real   sx, sy, s, t, cax, cay, s_len;

    sx = x2 - x1;
    sy = y2 - y1;
    cax = x - x1;
    cay = y - y1;

    s_len = sx * sx + sy * sy;
    if( s_len == 0.0 )
        return( FALSE );

    s = (cax * sx + cay * sy) / s_len;
    t = (cax * (-sy) + cay * sx) / s_len;

    *xwx1 = 1.0 - s;
    *xwy1 = t;
    *xwx2 = s;
    *xwy2 = -t;
    *ywx1 = -t;
    *ywy1 = 1.0 - s;
    *ywx2 = t;
    *ywy2 = s;

    return( TRUE );
}

public  BOOLEAN  get_prediction_weights_2d(
    Real   x,
    Real   y,
    int    n_points,
    Real   xs[],
    Real   ys[],
    Real   *x_weights[2],
    Real   *x_constant,
    Real   *y_weights[2],
    Real   *y_constant )
{
    int   dim, p, p1, p2, n_pairs;
    Real  xwx1, xwy1, xwx2, xwy2, ywx1, ywy1, ywx2, ywy2;

    *x_constant = 0.0;
    *y_constant = 0.0;

    for_less( dim, 0, 2 )
    {
        for_less( p, 0, n_points )
        {
            x_weights[dim][p] = 0.0;
            y_weights[dim][p] = 0.0;
        }
    }

    n_pairs = 0;
    for_less( p1, 0, n_points-1 )
    {
        for_less( p2, p1+1, n_points )
        {
            if( get_two_point_prediction( x, y, xs[p1], ys[p1], xs[p2], ys[p2],
                                          &xwx1, &xwy1, &xwx2, &xwy2,
                                          &ywx1, &ywy1, &ywx2, &ywy2 ) )
            {
                x_weights[0][p1] += xwx1;
                x_weights[1][p1] += xwy1;
                x_weights[0][p2] += xwx2;
                x_weights[1][p2] += xwy2;
                y_weights[0][p1] += ywx1;
                y_weights[1][p1] += ywy1;
                y_weights[0][p2] += ywx2;
                y_weights[1][p2] += ywy2;
                ++n_pairs;
            }
        }
    }

    for_less( dim, 0, 2 )
    {
        for_less( p, 0, n_points )
        {
            x_weights[dim][p] /= (Real) n_pairs;
            y_weights[dim][p] /= (Real) n_pairs;
        }
    }

#ifdef DEBUG
    test_solution( 0, x, y, n_points, xs, ys, x_weights[0], x_weights[1],
                   *x_constant );
    test_solution( 1, x, y, n_points, xs, ys, y_weights[0], y_weights[1],
                   *y_constant );
#endif

    return( TRUE );
}

#ifdef DEBUG
private  void  test_solution_3d(
    int    dim,
    Real   x,
    Real   y,
    Real   z,
    int    n_points,
    Real   xs[],
    Real   ys[],
    Real   zs[],
    Real   *weights[3] )
{
    int        iter, n_iters, p;
    Real       y_angle, z_angle, x_trans, y_trans, z_trans;
    Real       correct, value, ps[3];
    Transform  transform, y_rotation, z_rotation, translation;

    n_iters = 100;

    for_less( iter, 0, n_iters )
    {
        z_angle = 2.0 * PI * get_random_0_to_1();
        y_angle = 2.0 * PI * get_random_0_to_1();
        x_trans = 10.0 * get_random_0_to_1() - 5.0;
        y_trans = 10.0 * get_random_0_to_1() - 5.0;
        z_trans = 10.0 * get_random_0_to_1() - 5.0;

        make_rotation_transform( y_angle, Y, &y_rotation );
        make_rotation_transform( z_angle, Z, &z_rotation );
        make_translation_transform( x_trans, y_trans, z_trans, &translation );
        concat_transforms( &transform, &translation, &y_rotation );
        concat_transforms( &transform, &transform, &z_rotation );

        value = 0.0;

        for_less( p, 0, n_points )
        {
            transform_point( &transform, xs[p], ys[p], zs[p],
                             &ps[0], &ps[1], &ps[2] );
            value += weights[0][p] * ps[0];
            value += weights[1][p] * ps[1];
            value += weights[2][p] * ps[2];
        }

        transform_point( &transform, x, y, z, &ps[0], &ps[1], &ps[2] );
        correct = ps[dim];

        if( !numerically_close( value, correct, 1.0e-6 ) )
        {
            print( "get_prediction_weights_3d: %g %g\n", correct, value );
            break;
        }
    }
}
#endif

#define TOLERANCE  1.0e-1

private  BOOLEAN   get_four_point_prediction(
    Real   ax,
    Real   ay,
    Real   az,
    Real   ax1,
    Real   ay1,
    Real   az1,
    Real   ax2,
    Real   ay2,
    Real   az2,
    Real   ax3,
    Real   ay3,
    Real   az3,
    Real   ax4,
    Real   ay4,
    Real   az4,
    Real   weights[] )
{
    Real  x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;
    Real  x12, y12, z12, x23, x24, x34, x13, x14;
    Real  a, b, c, denom, v1_len, v2_len, max_len;

    x4 = ax - ax1;
    y4 = ay - ay1;
    z4 = az - az1;
    x1 = ax2 - ax1;
    y1 = ay2 - ay1;
    z1 = az2 - az1;
    x2 = ax3 - ax1;
    y2 = ay3 - ay1;
    z2 = az3 - az1;
    x3 = ax4 - ax1;
    y3 = ay4 - ay1;
    z3 = az4 - az1;

    v1_len = sqrt( x1 * x1 + y1 * y1 + z1 * z1 );
    v2_len = sqrt( x2 * x2 + y2 * y2 + z2 * z2 );

    max_len = MAX( v1_len, v2_len );

    x12 = x1 * y2 - x2 * y1;
    y12 = y1 * z2 - y2 * z1;
    z12 = z1 * x2 - z2 * x1;
    denom = x12 * z3 + z12 * y3 + y12 * x3;

    if( v1_len == 0.0 || v2_len == 0.0 ||
        FABS( denom / v1_len / v2_len / max_len ) < TOLERANCE )
        return( FALSE );

    x23 = x2 * y3 - x3 * y2;
    x24 = x2 * y4 - x4 * y2;
    x34 = x3 * y4 - x4 * y3;
    x13 = x1 * y3 - x3 * y1;
    x14 = x1 * y4 - x4 * y1;

    a = ( z4 * x23 + z3 * (-x24) + z2 * x34) / denom;
    b = (-z4 * x13 + z3 *   x14  - z1 * x34) / denom;
    c = ( z4 * x12 + z2 * (-x14) + z1 * x24) / denom;

    weights[0] = 1.0 - a - b - c;
    weights[1] = a;
    weights[2] = b;
    weights[3] = c;

    return( TRUE );
}

public  BOOLEAN  get_prediction_weights_3d(
    Real   x,
    Real   y,
    Real   z,
    int    n_points,
    Real   xs[],
    Real   ys[],
    Real   zs[],
    Real   *x_weights[3],
    Real   *y_weights[3],
    Real   *z_weights[3] )
{
    int   p, p1, p2, p3, p4, n_quads, dim;
    Real  weights4[4];

    for_less( p, 0, n_points )
    for_less( dim, 0, N_DIMENSIONS )
    {
        x_weights[dim][p] = 0.0;
        y_weights[dim][p] = 0.0;
        z_weights[dim][p] = 0.0;
    }

    n_quads = 0;
    for_less( p1, 0, n_points-3 )
    for_less( p2, p1+1, n_points-2 )
    for_less( p3, p2+1, n_points-1 )
    for_less( p4, p3+1, n_points )
    {
        if( get_four_point_prediction( x, y, z,
                                      xs[p1], ys[p1], zs[p1],
                                      xs[p2], ys[p2], zs[p2],
                                      xs[p3], ys[p3], zs[p3],
                                      xs[p4], ys[p4], zs[p4],
                                      weights4 ) )
        {
            x_weights[0][p1] += weights4[0];
            x_weights[0][p2] += weights4[1];
            x_weights[0][p3] += weights4[2];
            x_weights[0][p4] += weights4[3];
            y_weights[1][p1] += weights4[0];
            y_weights[1][p2] += weights4[1];
            y_weights[1][p3] += weights4[2];
            y_weights[1][p4] += weights4[3];
            z_weights[2][p1] += weights4[0];
            z_weights[2][p2] += weights4[1];
            z_weights[2][p3] += weights4[2];
            z_weights[2][p4] += weights4[3];
            ++n_quads;
        }
    }

    if( n_quads == 0 )
        return( FALSE );

    for_less( p, 0, n_points )
    for_less( dim, 0, N_DIMENSIONS )
    {
        x_weights[dim][p] /= (Real) n_quads;
        y_weights[dim][p] /= (Real) n_quads;
        z_weights[dim][p] /= (Real) n_quads;
    }

#ifdef DEBUG
    test_solution_3d( 0, x, y, z, n_points, xs, ys, zs, x_weights );
    test_solution_3d( 1, x, y, z, n_points, xs, ys, zs, y_weights );
    test_solution_3d( 2, x, y, z, n_points, xs, ys, zs, z_weights );
#endif

    return( TRUE );
}
