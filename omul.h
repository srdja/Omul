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
#ifndef __OMUL_H__
#define __OMUL_H__

#include <GL/gl.h>

void
omul_matrix_identity   (GLfloat matrix[]);

void
omul_matrix_multipy    (GLfloat out[],
                        GLfloat m1[],
                        GLfloat m2[]);

void
omul_matrix_multiply_ip(GLfloat m[],
                        GLfloat v[],
                        GLfloat off)

void
omul_rotate            (GLfloat matrix[],
                        GLfloat angle,
                        GLfloat x,
                        GLfloat y,
                        GLfloat z);

void
omul_multiply_vertex   (GLfloat out[],
                        GLfloat m[],
                        GLfloat v[],
                        GLint off);

void
omul_scale             (GLfloat matrix[],
                        GLfloat x,
                        GLfloat y,
                        GLfloat z);

void
omul_translate         (GLfloat matrix[],
                        GLfloat x,
                        GLfloat y,
                        GLfloat z);

void
omul_copy              (GLfloat cp[],
                        GLfloat matrix[]);

void
omul_look_at           (GLfloat matrix[],
                        GLfloat eye_x,
                        GLfloat eye_y,
                        GLfloat eye_z,
                        GLfloat center_x,
                        GLfloat center_y,
                        GLfloat center_z,
                        GLfloat up_x,
                        GLfloat up_y,
                        GLfloat up_z);

void
matrix_ortho           (GLfloat matrix[],
                        GLfloat left,
                        GLfloat right,
                        GLfloat bottom,
                        GLfloat top,
                        GLfloat n,
                        GLfloat f);

void
omul_normalize_vector  (GLfloat *x,
                        GLfloat *y,
                        GLfloat *z);

#endif
