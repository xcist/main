// Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include "TreePhantom.h"
#include "BaseObject.h"
#include "Cube.h"
#include "Sphere.h"
#include "Cylinder.h"
  
  extern "C" {
    int phantomlex(void);
    extern int linecount;
  }

  extern int phantomparse(void);
  extern int phantomdebug;
  
  // This contains the final stack of objects
  LinearPhantom ostack;
  BaseObject *top;
  int objectCount;
  MaterialTable localtab;
  
  // These are the state variables modified by the
  // assignments
  double lengthx, lengthy, lengthz;
  double radius1, radius2;
  double axis_x, axis_y, axis_z;
  double rho;
  double ax_x, ax_y, ax_z;
  double ay_x, ay_y, ay_z;
  double az_x, az_y, az_z;
  bool ax_given, ay_given, az_given, axis_given; 
  double x, y, z;
  double length, radius;
  double scale_factor; 
  int clip_count;
  double bb_radius; 
  double clipcoefs[50][4];
  std::string material;
  
  void phantomerror(const char *str)  {
    fprintf(stderr,"error: %s on line %d\n",str,linecount);
    exit(1);
  }
  
  int mapMaterialNameToOrdinal(std::string materialname) {
    MaterialTable::iterator result = std::find(localtab.begin(),
					       localtab.end(),
					       materialname);
    if (result == localtab.end()) {
      localtab.push_back(materialname);
      return localtab.size();
    }
    return result - localtab.begin() + 1;
  }

  void printMaterialTable(FILE *fp) {
    std::string myString;
    const char *ch;
    for(MaterialTable::iterator result = localtab.begin();result != localtab.end();result++)
      {
	myString = *result;
	ch = myString.data();
	fprintf(fp,"%s\n",ch);
      }
  }

  void resetParameters() {
    lengthx = lengthy = lengthz = 1;
    radius1 = radius2 = 1;
    axis_x = axis_y = axis_z = 0;
    rho = 1.0;
    ax_x = ax_y = ax_z = 0;
    ay_x = ay_y = ay_z = 0;
    az_x = az_y = az_z = 0;
    axis_given = false;
    ax_given = ay_given = az_given = false;
    x = y = z = 0;
    length = 1;
    radius = 1;
    top = NULL;
    clip_count = 0;
    bb_radius = 0;
    material = "water";
  }
  
  LinearPhantom parsePhantomDefFile(double scaleFactor, MaterialTable &dmtab) {
    ostack.clear();
    objectCount = 0;
    scale_factor = scaleFactor;
    resetParameters();
    localtab.clear();
    phantomparse();
    dmtab = localtab;
    return ostack;
  }

  void cross(double a1, double a2, double a3,
	     double b1, double b2, double b3,
	     double&c1, double&c2, double&c3) {
    c1 = a2*b3-a3*b2;
    c2 = a3*b1-a1*b3;
    c3 = a1*b2-a2*b1;
  }
  
  int sign(double x) {
    if (x >=0 ) return 1;
    if (x <=0 ) return -1;
    return 0;
  }
  
  void yyexpect(char *str) {
    fprintf(stderr,"Error - while parsing the phantom file, expected %s on line %d\n",str,linecount);
    exit(1);
  }
  
  %}

%union {
  double val;
  double vec[3];
  char str[128];
}

%token SPHERE BOX CYLX CYLY CYLZ CYL ELLIP ELLIPFREE
%token ELLIPCYL ECYLX ECYLY ECYLZ CONE CONEX CONEY CONEZ
%token RHO SIN COS SQRT
%token MATERIAL

%right '='
%left '-' '+'
%left '*' '/'
%left NEG

%token <val> LENGTHX LENGTHY LENGTHZ RADIUS1 NUMBER RADIUS2 
%token <vec> AXIS AX AY AZ 
%type  <vec> vectorExpression
%type  <val> expression
%token <str> STRING

%%
objectdefs:
| objectdefs objectdef
;

objectdef: 
'{' labeldef volumedef dmspec '}' 
{
  if (top != NULL) {
    top->Transform().Translate(Vec3(scale_factor*x,scale_factor*y,scale_factor*z));
    top->Density(rho);
    top->Priority(++objectCount);
    top->Material(mapMaterialNameToOrdinal(material));
    for (int m=0;m<clip_count;m++) 
      top->AddClipPlane(ClipPlane(Vec3(clipcoefs[m][0],clipcoefs[m][1],clipcoefs[m][2]),clipcoefs[m][3]));
    top->UpdateBoundingSphere();
    top->ID(ostack.size()+1);
    ostack.push_back(top);
  }
  resetParameters();
}
;


dmspec:
materialdef | densitydef | densitydef materialdef | materialdef densitydef
;

materialdef:
MATERIAL '=' STRING {material = std::string(strdup($3))} 
;

densitydef:
RHO '=' expression { rho = $3 }
;

labeldef:
| quotedstring
;

quotedstring:
'"' NUMBER '"'
;

volumedef: sphereDef
| boxDef
| cylinderXDef
| cylinderYDef
| cylinderZDef
| cylinderDef
| ellipsoidDef
| freeEllipsoidDef
| ellipticalCylDef
| eCylinderXDef
| eCylinderYDef
| eCylinderZDef
| coneDef
| coneXDef
| coneYDef
| coneZDef
;

commonAssignments: 
positionAssignment 
| clipPlanes;

sphereDef:
'[' SPHERE ':' sphereAssignmentBlock ']' {
  top = new Sphere;
  top->Transform().Scale(Vec3(scale_factor*radius,scale_factor*radius,scale_factor*radius));
}
;

sphereAssignmentBlock:
| sphereAssignmentBlock radiusAssignment
| sphereAssignmentBlock commonAssignments
;

boxDef:
'[' BOX ':' boxAssignmentBlock ']' {
  top = new Cube;
  top->Transform().Scale(Vec3(scale_factor*lengthx,scale_factor*lengthy,scale_factor*lengthz));
  if (!ax_given)
    cross(ay_x,ay_y,ay_z,az_x,az_y,az_z,ax_x,ax_y,ax_z);
  if (!ay_given)
    cross(az_x,az_y,az_z,ax_x,ax_y,ax_z,ay_x,ay_y,ay_z);
  if (!az_given)
    cross(ax_x,ax_y,ax_z,ay_x,ay_y,ay_z,az_x,az_y,az_z);
  double a1 = sqrt(ax_x*ax_x + ax_y*ax_y);
  double a2 = (az_x*ax_x + az_y*ax_y);
  double a3 = (az_y*ax_x - az_x*ax_y);
  double a4 = sqrt(a2*a2+az_z*az_z);
  double rotangle = 0;
  if (a3 != 0) {
    if (a4 == 0) {
      rotangle = 90.0*sign(a3);
    } else {
      rotangle = atan2(a3,a4*a1)*180.0/M_PI;
    }
  }
  top->Transform().RotateX(rotangle);
  rotangle = 0;
  if (ax_z != 0) {
    if (a1 == 0) {
      rotangle = 90*sign(ax_z);
    } else {
      rotangle = atan2(ax_z,a1)*180.0/M_PI;
    }
  }
  top->Transform().RotateY(rotangle);
  rotangle = 0;
  if (ax_y != 0) {
    if (ax_x == 0) {
      rotangle = 90*sign(ax_y);
    } else {
      rotangle = atan2(ax_y,ax_x)*180.0/M_PI;
    }
  }
  top->Transform().RotateZ(rotangle);
}
;

boxAssignmentBlock:
| boxAssignmentBlock lengthsAssignment
| boxAssignmentBlock axesAssignment
| boxAssignmentBlock commonAssignments
;

cylinderXDef:
'[' CYLX ':' acylinderAssignmentBlock ']' {
  top = new Cylinder;
  top->Transform().Scale(Vec3(scale_factor*radius,scale_factor*radius,scale_factor*length));
  top->Transform().RotateY(90);
}
;

cylinderYDef:
'[' CYLY ':' acylinderAssignmentBlock ']' {
  top = new Cylinder;
  top->Transform().Scale(Vec3(scale_factor*radius,scale_factor*radius,scale_factor*length));
  top->Transform().RotateX(90);
}
;

cylinderZDef:
'[' CYLZ ':' acylinderAssignmentBlock ']' {
  top = new Cylinder;
  top->Transform().Scale(Vec3(scale_factor*radius,scale_factor*radius,scale_factor*length));
  top->Transform().RotateZ(90);
}
;

acylinderAssignmentBlock:
| acylinderAssignmentBlock radiusAssignment 
| acylinderAssignmentBlock lengthAssignment
| acylinderAssignmentBlock commonAssignments
;

cylinderDef:
'[' CYL ':' cylinderAssignmentBlock ']' {
  top = new Cylinder;
  top->Transform().Scale(Vec3(scale_factor*radius,scale_factor*radius,scale_factor*length));
   double rotang;
  rotang = -atan2(sqrt(axis_x*axis_x+axis_y*axis_y),axis_z)*180.0/M_PI;
  top->Transform().RotateX(rotang);
  rotang = -atan2(axis_x,axis_y)*180.0/M_PI;
  top->Transform().RotateZ(rotang);
}
;

cylinderAssignmentBlock:
| cylinderAssignmentBlock radiusAssignment 
| cylinderAssignmentBlock lengthAssignment
| cylinderAssignmentBlock axisAssignment
| cylinderAssignmentBlock commonAssignments
;

ellipsoidDef:
'[' ELLIP ':' ellipsoidAssignmentBlock ']' {
  top = new Sphere;
  top->Transform().Scale(Vec3(scale_factor*lengthx,scale_factor*lengthy,scale_factor*lengthz));
}
;

ellipsoidAssignmentBlock:
| ellipsoidAssignmentBlock ellipsoidAssignmentStatement;

ellipsoidAssignmentStatement:
lengthsAssignment 
| commonAssignments
;

freeEllipsoidDef:
'[' ELLIPFREE ':' freeEllipsoidAssignmentBlock ']' { 
  top = new Sphere;
  top->Transform().Scale(Vec3(scale_factor*lengthx,scale_factor*lengthy,scale_factor*lengthz));
  if (!ax_given)
    cross(ay_x,ay_y,ay_z,az_x,az_y,az_z,ax_x,ax_y,ax_z);
  if (!ay_given)
    cross(az_x,az_y,az_z,ax_x,ax_y,ax_z,ay_x,ay_y,ay_z);
  if (!az_given)
    cross(ax_x,ax_y,ax_z,ay_x,ay_y,ay_z,az_x,az_y,az_z);
  double a1 = sqrt(ax_x*ax_x + ax_y*ax_y);
  double a2 = (az_x*ax_x + az_y*ax_y);
  double a3 = (az_y*ax_x - az_x*ax_y);
  double a4 = sqrt(a2*a2+az_z*az_z);
  double rotangle = 0;
  if (a3 != 0) {
    if (a4 == 0) {
      rotangle = 90.0*sign(a3);
    } else {
      rotangle = atan2(a3,a4*a1)*180.0/M_PI;
    }
  }
  top->Transform().RotateX(rotangle);
  rotangle = 0;
  if (ax_z != 0) {
    if (a1 == 0) {
      rotangle = 90*sign(ax_z);
    } else {
      rotangle = atan2(ax_z,a1)*180.0/M_PI;
    }
  }
  top->Transform().RotateY(rotangle);
  rotangle = 0;
  if (ax_y != 0) {
    if (ax_x == 0) {
      rotangle = 90*sign(ax_y);
    } else {
      rotangle = atan2(ax_y,ax_x)*180.0/M_PI;
    }
  }
  top->Transform().RotateZ(rotangle);
};

freeEllipsoidAssignmentBlock:
| freeEllipsoidAssignmentBlock lengthsAssignment
| freeEllipsoidAssignmentBlock axesAssignment
| freeEllipsoidAssignmentBlock commonAssignments
;


ellipticalCylDef:
'[' ELLIPCYL ':' eCylinderAssignmentBlock ']' {
  top = new Cylinder;
  top->Transform().Scale(Vec3(scale_factor*lengthx,scale_factor*lengthy,scale_factor*length));
  if (axis_given) {
    az_x = axis_x;
    az_y = axis_y;
    az_z = axis_z;
  }
  if (!ax_given)
    cross(ay_x,ay_y,ay_z,az_x,az_y,az_z,ax_x,ax_y,ax_z);
  if (!ay_given)
    cross(az_x,az_y,az_z,ax_x,ax_y,ax_z,ay_x,ay_y,ay_z);
  if (!az_given)
    cross(ax_x,ax_y,ax_z,ay_x,ay_y,ay_z,az_x,az_y,az_z);
  double a1 = sqrt(ax_x*ax_x + ax_y*ax_y);
  double a2 = (az_x*ax_x + az_y*ax_y);
  double a3 = (az_y*ax_x - az_x*ax_y);
  double a4 = sqrt(a2*a2+az_z*az_z);
  double rotangle = 0;
  if (a3 != 0) {
    if (a4 == 0) {
      rotangle = 90.0*sign(a3);
    } else {
      rotangle = atan2(a3,a4*a1)*180.0/M_PI;
    }
  }
  top->Transform().RotateX(rotangle);
  rotangle = 0;
  if (ax_z != 0) {
    if (a1 == 0) {
      rotangle = 90*sign(ax_z);
    } else {
      rotangle = atan2(ax_z,a1)*180.0/M_PI;
    }
  }
  top->Transform().RotateY(rotangle);
  rotangle = 0;
  if (ax_y != 0) {
    if (ax_x == 0) {
      rotangle = 90*sign(ax_y);
    } else {
      rotangle = atan2(ax_y,ax_x)*180.0/M_PI;
    }
  }
  top->Transform().RotateZ(rotangle);
};

eCylinderAssignmentBlock:
| eCylinderAssignmentBlock axisAssignment
| eCylinderAssignmentBlock axesAssignment
| eCylinderAssignmentBlock lengthsAssignment
| eCylinderAssignmentBlock lengthAssignment
| eCylinderAssignmentBlock commonAssignments
;

eCylinderXDef:
'[' ECYLX ':' aECylAssignmentBlock ']' {
  top = new Cylinder();
  top->Transform().Scale(Vec3(scale_factor*lengthz,scale_factor*lengthy,scale_factor*length));
  top->Transform().RotateY(90);
}
;

eCylinderYDef:
'[' ECYLY ':' aECylAssignmentBlock ']'  {
  top = new Cylinder();
  top->Transform().Scale(Vec3(scale_factor*lengthx,scale_factor*lengthz,scale_factor*length));
  top->Transform().RotateX(90);
}
;

eCylinderZDef:
'[' ECYLZ ':' aECylAssignmentBlock ']'  {
  top = new Cylinder();
  top->Transform().Scale(Vec3(scale_factor*lengthx,scale_factor*lengthy,scale_factor*length));
}
;

aECylAssignmentBlock:
| aECylAssignmentBlock lengthAssignment
| aECylAssignmentBlock lengthsAssignment
| aECylAssignmentBlock commonAssignments
;

coneDef:
'[' CONE ':' coneAssignmentBlock ']'
;

coneAssignmentBlock:
| coneAssignmentBlock axisAssignment
| coneAssignmentBlock lengthAssignment
| coneAssignmentBlock radiiAssignment
| coneAssignmentBlock commonAssignments
;

coneXDef:
'[' CONEX ':' aConeAssignmentBlock ']'
;

coneYDef:
'[' CONEY ':' aConeAssignmentBlock ']'
;

coneZDef:
'[' CONEZ ':' aConeAssignmentBlock ']'
;

aConeAssignmentBlock:
| coneAssignmentBlock lengthAssignment
| coneAssignmentBlock radiiAssignment
| coneAssignmentBlock commonAssignments
;


// Position the center of the object
positionAssignment:
'x' '=' expression   { x = $3 }
| 'y' '=' expression { y = $3 }
| 'z' '=' expression { z = $3 }
;

// Set the radius of the object
radiusAssignment:
'r' '=' expression { radius = $3 }
;

// Set the lengths of an object along the x, y & z axes
lengthsAssignment:
LENGTHX '=' expression   { lengthx = $3 }
| LENGTHY '=' expression { lengthy = $3 }
| LENGTHZ '=' expression { lengthz = $3 }
;

// Set a single length
lengthAssignment:
'l' '=' expression     { length = $3 }  
;

// Set multiple radii
radiiAssignment:
RADIUS1 '=' expression   { radius1 = $3 }
| RADIUS2 '=' expression { radius2 = $3 }
;

axisAssignment:
AXIS vectorExpression  { axis_x = $2[0]; axis_y = $2[1]; axis_z = $2[2]; axis_given = true; }
;

axesAssignment:
AX vectorExpression  { ax_x = $2[0]; ax_y = $2[1]; ax_z = $2[2]; ax_given = true; }
| AY vectorExpression  { ay_x = $2[0]; ay_y = $2[1]; ay_z = $2[2]; ay_given = true;}
| AZ vectorExpression  { az_x = $2[0]; az_y = $2[1]; az_z = $2[2]; az_given = true;}
;

vectorExpression:
'(' expression ',' expression ',' expression ')'
{ $$[0] = $2; $$[1] = $4; $$[2] = $6; }
;

clipPlanes:
//	  clipPlanes clipPlane
//	  ;
//
//clipPlane:
'x' '<' expression {
  clipcoefs[clip_count][0] = -1.0;
  clipcoefs[clip_count][1] = 0.0;
  clipcoefs[clip_count][2] = 0.0;
  clipcoefs[clip_count][3] = -scale_factor*$3;
  clip_count++;
}
| 'x' '>' expression {
  clipcoefs[clip_count][0] = 1.0;
  clipcoefs[clip_count][1] = 0.0;
  clipcoefs[clip_count][2] = 0.0;
  clipcoefs[clip_count][3] = scale_factor*$3;
  clip_count++;
}
| 'y' '<' expression {
  clipcoefs[clip_count][0] = 0.0;
  clipcoefs[clip_count][1] = -1.0;
  clipcoefs[clip_count][2] = 0.0;
  clipcoefs[clip_count][3] = -scale_factor*$3;
  clip_count++;
}
| 'y' '>' expression {
  clipcoefs[clip_count][0] = 0.0;
  clipcoefs[clip_count][1] = 1.0;
  clipcoefs[clip_count][2] = 0.0;
  clipcoefs[clip_count][3] = scale_factor*$3;
  clip_count++;
}
| 'z' '<' expression {
  clipcoefs[clip_count][0] = 0.0;
  clipcoefs[clip_count][1] = 0.0;
  clipcoefs[clip_count][2] = -1.0;
  clipcoefs[clip_count][3] = -scale_factor*$3;
  clip_count++;
}
| 'z' '>' expression {
  clipcoefs[clip_count][0] = 0.0;
  clipcoefs[clip_count][1] = 0.0;
  clipcoefs[clip_count][2] = 1.0;
  clipcoefs[clip_count][3] = scale_factor*$3;
  clip_count++;
}
| 'r' vectorExpression '<' expression {
  double norm = sqrt($2[0]*$2[0] + $2[1]*$2[1] + $2[2]*$2[2]);
  clipcoefs[clip_count][0] = -$2[0]/norm;
  clipcoefs[clip_count][1] = -$2[1]/norm;
  clipcoefs[clip_count][2] = -$2[2]/norm;
  clipcoefs[clip_count][3] = -scale_factor*$4;
  clip_count++;
}
| 'r' vectorExpression '>' expression {
  double norm = sqrt($2[0]*$2[0] + $2[1]*$2[1] + $2[2]*$2[2]);
  clipcoefs[clip_count][0] = $2[0]/norm;
  clipcoefs[clip_count][1] = $2[1]/norm;
  clipcoefs[clip_count][2] = $2[2]/norm;
  clipcoefs[clip_count][3] = scale_factor*$4;
  clip_count++;
}
;

expression:
NUMBER
| SIN '(' expression ')' { $$ = sin($3*M_PI/180.0) }
| COS '(' expression ')' { $$ = cos($3*M_PI/180.0) }
| SQRT '(' expression ')' { $$ = sqrt($3) }
| expression '+' expression { $$ = $1 + $3 }
| expression '-' expression { $$ = $1 - $3 }
| expression '*' expression { $$ = $1 * $3 }
| expression '/' expression { $$ = $1 / $3 }
| '-' expression %prec NEG  { $$ = -$2 }
| '(' expression ')'        { $$ = $2 }
;
