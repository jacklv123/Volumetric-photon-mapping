#include "object.h"

CObj::CObj(int tp, int idx, float *cod) {

  index = idx;
  optic = OPT_NONE;
  refractive = 1.0;
  type = tp;

  int ncods = 0;
  switch(tp) {
  case TYPE_SPHERE : ncods=4; break;
  case TYPE_PLANE  : ncods=2; break;
  }
  for(int i=0; i<ncods; i++) coords[i] = cod[i];
  for(int i=0; i<3; i++) color[i] = 1.0;
}


Vector3
CObj::calcSphereNormal(const Vector3 &P, const Vector3 &O)
{
  //球の中心を引き算-->法線ベクトル
  Vector3 center(coords[0], coords[1], coords[2]);
  Vector3 ans = P - center;
  ans.normalize();
  return ans;
}

Vector3
CObj::calcPlaneNormal(const Vector3 &P, const Vector3 &O)
{
  int axis = (int) this->coords[0];
  float N[3] = {0.0,0.0,0.0};
  N[axis] = O[axis] - this->coords[1];      //Vector From Surface to Light
  Vector3 ans(N);
  ans.normalize();
  return ans;
}


double
CObj::calcSphereIntersection(const Vector3 &r, const Vector3 &o) //Ray-Sphere Intersection: r=Ray Direction, o=Ray Origin
{
  //  s = Sphere Center Translated into Coordinate Frame of Ray Origin
  Vector3 center(coords[0], coords[1], coords[2]);
  Vector3 s = center - o;

  //radius=Sphere Radius
  float radius = coords[3];

  //Intersection of Sphere and Line     =       Quadratic Function of Distance
  float A = dot(r,r);
  float B = -2.0 * dot(s,r);
  float C = dot(s,s) - radius * radius;
  float D = B * B - 4 * A * C;

  //  二次方程式に解がある場合のみ，距離が算出できる
  if (D > 0.0) {
    float sign = (C < -0.00001) ? 1 : -1;
    return (-B + sign*sqrt(D))/(2*A);
  }
  return NOT_INTERSECTED;
}

double
CObj::calcPlaneIntersection(const Vector3 &r, const Vector3 &o)
{
  int axis = (int) coords[0];            //Determine Orientation of Axis-Aligned Plane
  if (r[axis] != 0.0){                        //Parallel Ray -> No Intersection
    return  (coords[1] - o[axis]) / r[axis]; //Solve Linear Equation (rx = p-o)
  }
  return NOT_INTERSECTED;
}
