
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
