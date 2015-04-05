/* ----------------------------------------------------------------------------
 * Omul
 * Copyright (c) 2015 Srđan Panić
 *
 * This file is a part of OMUL
 *
 * This software is provided 'as-is', without any express or implied
 * warranty. In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgement in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 *
 * --------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <GL/gl.h>

#include "omul.h"

#define BLOCK_SIZE sizeof(GLfloat) * 16

/*-----------------------------------------------------------------------------
 * Makes the matrix m an identity matrix.
 *
 * @param[in, out] the matrix on which this operation is performed.
 *---------------------------------------------------------------------------*/
void
omul_matrix_identity(GLfloat m[])
{
    memset(m, 0, BLOCK_SIZE);

    m[0 ] = 1.0f;
    m[5 ] = 1.0f;
    m[10] = 1.0f;
    m[15] = 1.0f;
}

/*----------------------------------------------------------------------------
 * Multiplies two 4x4 matrices and stores product in the third matrix array.
 *
 * @param[out] out the output matrix in which the result of the multiplication is
 *                 stored
 * @param[in] m1 matrix that is being multiplied by m2
 * @param[in] m2 matrix that multiplies m1
 *---------------------------------------------------------------------------*/
void
omul_matrix_multiply(GLfloat out[], GLfloat m1[], GLfloat m2[])
{
    out[0 ] = m1[0 ] * m2[0 ] + m1[4 ] * m2[1 ] + m1[8 ] * m2[2 ] + m1[12] * m2[3];
    out[1 ] = m1[1 ] * m2[0 ] + m1[5 ] * m2[1 ] + m1[9 ] * m2[2 ] + m1[13] * m2[3];
    out[2 ] = m1[2 ] * m2[0 ] + m1[6 ] * m2[1 ] + m1[10] * m2[2 ] + m1[14] * m2[3];
    out[3 ] = m1[3 ] * m2[0 ] + m1[7 ] * m2[1 ] + m1[11] * m2[2 ] + m1[15] * m2[3];

    out[4 ] = m1[0 ] * m2[4 ] + m1[4 ] * m2[5 ] + m1[8 ] * m2[6 ] + m1[12] * m2[7];
    out[5 ] = m1[1 ] * m2[4 ] + m1[5 ] * m2[5 ] + m1[9 ] * m2[6 ] + m1[13] * m2[7];
    out[6 ] = m1[2 ] * m2[4 ] + m1[6 ] * m2[5 ] + m1[10] * m2[6 ] + m1[14] * m2[7];
    out[7 ] = m1[3 ] * m2[4 ] + m1[7 ] * m2[5 ] + m1[11] * m2[6 ] + m1[15] * m2[7];

    out[8 ] = m1[0 ] * m2[8 ] + m1[4 ] * m2[9 ] + m1[8 ] * m2[10] + m1[12] * m2[11];
    out[9 ] = m1[1 ] * m2[8 ] + m1[5 ] * m2[9 ] + m1[9 ] * m2[10] + m1[13] * m2[11];
    out[10] = m1[2 ] * m2[8 ] + m1[6 ] * m2[9 ] + m1[10] * m2[10] + m1[14] * m2[11];
    out[11] = m1[3 ] * m2[8 ] + m1[7 ] * m2[9 ] + m1[11] * m2[10] + m1[15] * m2[11];

    out[12] = m1[0 ] * m2[12] + m1[4 ] * m2[13] + m1[8 ] * m2[14] + m1[12] * m2[15];
    out[13] = m1[1 ] * m2[12] + m1[5 ] * m2[13] + m1[9 ] * m2[14] + m1[13] * m2[15];
    out[14] = m1[2 ] * m2[12] + m1[6 ] * m2[13] + m1[10] * m2[14] + m1[14] * m2[15];
    out[15] = m1[3 ] * m2[12] + m1[7 ] * m2[13] + m1[11] * m2[14] + m1[15] * m2[15];
}

/*-----------------------------------------------------------------------------
 * Multiplies a row major vertex with a column major matrix
 *
 * @param[out] out a 3 element vector in which the transformed vertex is stored
 * @param[in]  m the matrix by which the vertex is multiplied
 * @param[in]  vert the vertex being transformed
 * @param[in]  vert_offset the offset into the vertex array
 *---------------------------------------------------------------------------*/
void
omul_multiply_vertex(GLfloat out[], GLfloat m[], GLfloat v[], GLint off)
{
    /* implicit w = 1 */
    out[0] = m[0] * v[off] + m[4] * v[off + 1] + m[8 ] * v[off + 2] + m[12];
    out[1] = m[1] * v[off] + m[5] * v[off + 1] + m[9 ] * v[off + 2] + m[13];
    out[2] = m[2] * v[off] + m[6] * v[off + 1] + m[10] * v[off + 2] + m[14];
}


/*-----------------------------------------------------------------------------
 * "In place" (not really) matrix-vertex multiplication.
 *
 * @param[in]  m the matrix that multiplies the vertex
 * @param[out] v the vertex being multiplied
 * @param[in]  off offset into the vertex array
 *---------------------------------------------------------------------------*/
void
omul_multiply_vertex_ip(GLfloat m[], GLfloat v[], GLint off)
{
    /* implicit w = 1 */
    float x = v[off    ];
    float y = v[off + 1];
    float z = v[off + 2];

    v[off    ] = m[0] * x + m[4] * y + m[8 ] * z + m[12];
    v[off + 1] = m[1] * x + m[5] * y + m[9 ] * z + m[13];
    v[off + 2] = m[2] * x + m[6] * y + m[10] * z + m[14];
}

/*-----------------------------------------------------------------------------
 * Multiplies the matrix m by the scale factor xyz
 *
 * @param[in, out] the matrix that is being multiplied
 * @param[in] x,y,z  the scale factors
 *---------------------------------------------------------------------------*/
void
omul_scale(GLfloat m[], GLfloat x, GLfloat y, GLfloat z)
{
    m[0 ] *= x;
    m[1 ] *= x;
    m[2 ] *= x;
    m[3 ] *= x;

    m[4 ] *= y;
    m[5 ] *= y;
    m[6 ] *= y;
    m[7 ] *= y;

    m[8 ] *= z;
    m[9 ] *= z;
    m[10] *= z;
    m[11] *= z;
}

/*-----------------------------------------------------------------------------
 * Multiplies the matrix m by the translation vector xyz.
 *
 * @param[in] m the matrix that is being multiplied
 * @param[in] (x, y, z) the translation vector
 *---------------------------------------------------------------------------*/
void
omul_translate(GLfloat m[], GLfloat x, GLfloat y, GLfloat z)
{
    m[12] += m[0] * x + m[4] * y + m[8 ] * z;
    m[13] += m[1] * x + m[5] * y + m[9 ] * z;
    m[14] += m[2] * x + m[6] * y + m[10] * z;
    m[15] += m[3] * x + m[7] * y + m[11] * z;
}

/*-----------------------------------------------------------------------------
 * Based on: https://www.opengl.org/sdk/docs/man2/xhtml/glRotate.xml spec.
 *
 * @param[in, out]  m the matrix that is being transformed
 * @param[in] angle rotation angle in degrees
 * @param[in] x, y, z Represent the axis vector around which the rotation is
 *            performed.
 *---------------------------------------------------------------------------*/
void
omul_rotate(GLfloat m[], GLfloat angle, GLfloat x, GLfloat y, GLfloat z)
{
    GLfloat rad   = (GLfloat) (angle * M_PI / 180.0f);
    GLfloat a_sin = (GLfloat) sin(rad);
    GLfloat a_cos = (GLfloat) cos(rad);

    omul_normalize_vector(&x , &y, &z);

    GLfloat xx = x * x;
    GLfloat xy = x * y;
    GLfloat xz = x * z;
    GLfloat yy = y * y;
    GLfloat yz = y * z;
    GLfloat zz = z * z;

    GLfloat c  = (GLfloat) (1.0f - a_cos);

    m[0 ] = xx * c + a_cos;
    m[1 ] = xy * c - (z * a_sin);
    m[2 ] = xz * c + (y * a_sin);

    m[4 ] = xy * c + (z * a_sin);
    m[5 ] = yy * c + a_cos;
    m[6 ] = yz * c - (x * a_sin);

    m[8 ] = xz * c - (y * a_sin);
    m[9 ] = yz * c + (x * a_sin);
    m[10] = zz * c + a_cos;
}

/*-----------------------------------------------------------------------------
 * Based on: https://www.opengl.org/sdk/docs/man2/xhtml/gluLookAt.xml
 *
 * @param[in, out] m
 * @param[in] eye_x x eye point
 * @param[in] eye_y y eye point
 * @param[in] eye_z z eye point
 * @param[in] center_x x reference point
 * @param[in] center_y y reference point
 * @param[in] center_z z reference point
 * @param[in] up_x direction vector x
 * @param[in] up_y direction vector y
 * @param[in] up_z direction vector z
 *---------------------------------------------------------------------------*/
void
omul_look_at(GLfloat m[],
             GLfloat eye_x,
             GLfloat eye_y,
             GLfloat eye_z,
             GLfloat center_x,
             GLfloat center_y,
             GLfloat center_z,
             GLfloat up_x,
             GLfloat up_y,
             GLfloat up_z)
{
    GLfloat forward_x = center_x - eye_x;
    GLfloat forward_y = center_y - eye_y;
    GLfloat forward_z = center_z - eye_z;

    GLfloat norm_fw = 1.0f / sqrt(forward_x * forward_x +
                                  forward_y * forward_y +
                                  forward_z * forward_z);

    forward_x *= norm_fw;
    forward_y *= norm_fw;
    forward_z *= norm_fw;

    GLfloat side_x_cp = forward_y * up_z - forward_z * up_y;
    GLfloat side_y_cp = forward_z * up_x - forward_x * up_z;
    GLfloat side_z_cp = forward_x * up_y - forward_y * up_x;

    GLfloat norm_s = 1.0f / sqrt(side_x_cp * side_x_cp +
                                 side_y_cp * side_y_cp +
                                 side_z_cp * side_z_cp);

    side_x_cp *= norm_s;
    side_y_cp *= norm_s;
    side_z_cp *= norm_s;

    GLfloat ux = side_y_cp * forward_z - side_z_cp * forward_y;
    GLfloat uy = side_z_cp * forward_x - side_x_cp * forward_z;
    GLfloat uz = side_x_cp * forward_y - side_y_cp * forward_x;

    m[0 ] = side_x_cp;
    m[1 ] = ux;
    m[2 ] = -forward_x;
    m[3 ] = 0;

    m[4 ] = side_y_cp;
    m[5 ] = uy;
    m[6 ] = -forward_y;
    m[7 ] = 0;

    m[8 ] = side_z_cp;
    m[9 ] = uz;
    m[10] = -forward_z;
    m[11] = 0;

    m[12] = 0;
    m[13] = 0;
    m[14] = 0;
    m[15] = 1.0f;

    matrix_translate(m, -eye_x, -eye_y, -eye_z);
}

/*-----------------------------------------------------------------------------
 * Based on: https://www.opengl.org/sdk/docs/man2/xhtml/glOrtho.xml
 *
 * @param[in, out] m the matrix
 * @param[in] left   coordinate for the left vertical clipping plane
 * @param[in] right  coordinate for the right vertical clipping plane
 * @param[in] bottom coordinate for the bottom horizontal clipping plane
 * @param[in] top    coordinate for the top horizontal clipping plane
 * @param[in] n      distance to the near depth clipping plane
 * @param[in] f      distance to the far depth clipping plane
 *---------------------------------------------------------------------------*/
void
omul_ortho(GLfloat m[],
           GLfloat left,
           GLfloat right,
           GLfloat bottom,
           GLfloat top,
           GLfloat n,
           GLfloat f)
{
    GLfloat x =  2.0f / (right - left);
    GLfloat y =  2.0f / (top - bottom);
    GLfloat z = -2.0f / (f - n);

    GLfloat tx = -((right + left) / (right - left));
    GLfloat ty = -((top + bottom) / (top - bottom));
    GLfloat tz = -((f + n)  / (f - n));

    memset(m, 0, BLOCK_SIZE);

    m[0 ] = x;
    m[5 ] = y;
    m[10] = z;

    m[12] = tx;
    m[13] = ty;
    m[14] = tz;
    m[15] = 1.0f;
}

/*-----------------------------------------------------------------------------
 * Copies a 4x4 matrix m into matrix cp
 *
 * @param[out] cp the copy of the matrix m
 * @param[in]  m the matrix being copied
 *---------------------------------------------------------------------------*/
inline void
omul_copy(GLfloat cp[], GLfloat m[])
{
    memcpy(cp, m, BLOCK_SIZE);
}

/*-----------------------------------------------------------------------------
 * Normalizes the vector x,y,z
 *
 * @param[in, out] x
 * @param[in, out] y
 * @param[in, out] z
 *---------------------------------------------------------------------------*/
inline void
omul_normalize_vector(GLfloat *x, GLfloat *y, GLfloat *z)
{
    GLfloat l = sqrt((*x) * (*x) +
                     (*y) * (*y) +
                     (*z) * (*z));

    if (l == 1)
        return;

    *x *= 1.0f / l;
    *y *= 1.0f / l;
    *z *= 1.0f / l;
}
