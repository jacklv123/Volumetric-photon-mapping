//CObj.h
#include <cstdlib>
#include "vector3.h"

using WebCore::Vector3;

#define TYPE_SPHERE   0
#define TYPE_PLANE    1
#define TYPE_TRIANGLE 2
#define TYPE_VOLUME 3

#define OPT_NONE 0
#define OPT_REFLECT 1
#define OPT_REFRACT 2

#define NOT_INTERSECTED (1.0e6)

class CObj {
  /**ŠÖ”**/
  public :
    CObj(int type, int idx, float *cod);
    int    getType()   { return type; }
    int    getOptics() { return optic; }
    int    getIndex()  { return index; }
    float  getRefractive()  { return refractive; }
    void   setOptics(int op)   { optic = op; }
    void   setRefractive(float rf)   { this->refractive = rf; }
    void   setColor(const float *cl) {
      for(int i=0; i<3; i++) color[i] = cl[i];
    }
    float coords[9];
    float color[3];

    Vector3 calcSphereNormal(const Vector3 &P, const Vector3 &O);
    Vector3 calcPlaneNormal(const Vector3 &P, const Vector3 &O);

    double  calcSphereIntersection(const Vector3 &ray, const Vector3 &org);
    double  calcPlaneIntersection(const Vector3 &ray, const Vector3 &org);

    /**•Ï”**/
  private :
    int type;
    int optic;
    int index;
    float refractive;
};

typedef struct SIntersectionStat {
  CObj    *obj;
  double  dist;
  SIntersectionStat() {
    dist = NOT_INTERSECTED;
    obj = NULL;
  }
} SIntersectionStat;
