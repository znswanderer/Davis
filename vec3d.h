#ifndef VEC_H
#define VEC_H

typedef struct {
  double x;
  double y;
  double z;
} vec;

#define VEC_ZERO (vec){ 0.0, 0.0, 0.0 }; 

#define vec_inc(v, w) { \
  (v).x += (w).x;	\
  (v).y += (w).y;	\
  (v).z += (w).z;	\
  }

#define vec_dec(v, w) { \
    (v).x -= (w).x;	\
    (v).y -= (w).y;	\
    (v).z -= (w).z;	\
  }

#define vec_add(v, w) (vec){ (v).x+(w).x, (v).y+(w).y, (v).z+(w).z }
#define vec_sub(v, w) (vec){ (v).x-(w).x, (v).y-(w).y, (v).z-(w).z }
#define vec_sclr(k, v) (vec){ (k)*(v).x, (k)*(v).y, (k)*(v).z }
#define vec_map(f, v) (vec){ f((v).x), f((v).y), f((v).z) }
#define vec_magnitude2(v) ((v).x*(v).x + (v).y*(v).y + (v).z*(v).z)
#define vec_magnitude(v) (sqrt(vec_magnitude2((v))))
#define vec_dot(v, w) ((v).x*(w).x + (v).y*(w).y + (v).z*(w).z)

#define vec_cross(a, b) (vec){			\
    a.y*b.z - a.z*b.y,	 a.z*b.x - a.x*b.z,   a.x*b.y - a.y*b.x }			


    



#endif
