%{
/* -----------------------------------------------------------------------
*   Program Name: phantom.l            
*   Authors:  Samit Basu, Bruno De Man, Jed Pack (GE Global Research)
*   Organization: GE Global Research and GE Healthcare, General Electric Company
*   Version:  6.0.3
*   Date:  Feb 3, 2015
* -----------------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include "phantom.tab.hpp"
// added by wrr/ejn
int linecount = 0;
extern YYSTYPE phantomlval;
#include "stdlib.h"
#include "math.h"
%}

INT				[0-9]+

%%
Sphere				return SPHERE;
Box				return BOX;
Cylinder_x			return CYLX;
Cylinder_y			return CYLY;
Cylinder_z			return CYLZ;
Cylinder			return CYL;
Ellipsoid			return ELLIP;
Ellipsoid_free			return ELLIPFREE;
Ellipt_Cyl			return ELLIPCYL;
Ellipt_Cyl_x			return ECYLX;
Ellipt_Cyl_y			return ECYLY;
Ellipt_Cyl_z			return ECYLZ;
Cone				return CONE;
Cone_x				return CONEX;
Cone_y				return CONEY;
Cone_z				return CONEZ;
dx				return LENGTHX;
dy				return LENGTHY;
dz				return LENGTHZ;
r1				return RADIUS1;
r2				return RADIUS2;
axis				return AXIS;
achse				return AXIS;
a_x				return AX;
a_y				return AY;
a_z				return AZ;
dichte				return RHO;
density				return RHO;
rho				return RHO;
material                        return MATERIAL;
{INT}("."{INT}?)? 		phantomlval.val = atof(phantomtext); return NUMBER;
"."{INT}			phantomlval.val = atof(phantomtext); return NUMBER;
[\'][^\'\n\r]*[\']              {
strcpy(phantomlval.str,phantomtext+1); 
phantomlval.str[strlen(phantomlval.str)-1] = 0;
return STRING;
}
[ \t]+				/* Skip */;
sin				return SIN;
cos				return COS;
sqrt				return SQRT;
#[^\n]*				/* Skip */;
\n				linecount++; 
.				return (int) phantomtext[0];

%%

int phantomwrap()
{
  return 1;
}
