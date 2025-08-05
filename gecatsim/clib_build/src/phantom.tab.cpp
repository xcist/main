/* -----------------------------------------------------------------------
*   Program Name: phantom.tab.cpp           
*   Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license
* -----------------------------------------------------------------------*/
/* A Bison parser, made by GNU Bison 1.875c.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* Written by Richard Stallman by simplifying the original so called
   ``semantic'' parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0

/* If NAME_PREFIX is specified substitute the variables and functions
   names.  */
#define yyparse phantomparse
#define yylex   phantomlex
#define yyerror phantomerror
#define yylval  phantomlval
#define yychar  phantomchar
#define yydebug phantomdebug
#define yynerrs phantomnerrs

#ifdef WIN32
#define _USE_MATH_DEFINES
typedef struct phantomlval {
  double val;
  double vec[3];
  char str[128];
} PHANTOMLVAL;

#endif

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     SPHERE = 258,
     BOX = 259,
     CYLX = 260,
     CYLY = 261,
     CYLZ = 262,
     CYL = 263,
     ELLIP = 264,
     ELLIPFREE = 265,
     ELLIPCYL = 266,
     ECYLX = 267,
     ECYLY = 268,
     ECYLZ = 269,
     CONE = 270,
     CONEX = 271,
     CONEY = 272,
     CONEZ = 273,
     RHO = 274,
     SIN = 275,
     COS = 276,
     SQRT = 277,
     MATERIAL = 278,
     NEG = 279,
     LENGTHX = 280,
     LENGTHY = 281,
     LENGTHZ = 282,
     RADIUS1 = 283,
     NUMBER = 284,
     RADIUS2 = 285,
     AXIS = 286,
     AX = 287,
     AY = 288,
     AZ = 289,
     STRING = 290
   };
#endif
#define SPHERE 258
#define BOX 259
#define CYLX 260
#define CYLY 261
#define CYLZ 262
#define CYL 263
#define ELLIP 264
#define ELLIPFREE 265
#define ELLIPCYL 266
#define ECYLX 267
#define ECYLY 268
#define ECYLZ 269
#define CONE 270
#define CONEX 271
#define CONEY 272
#define CONEZ 273
#define RHO 274
#define SIN 275
#define COS 276
#define SQRT 277
#define MATERIAL 278
#define NEG 279
#define LENGTHX 280
#define LENGTHY 281
#define LENGTHZ 282
#define RADIUS1 283
#define NUMBER 284
#define RADIUS2 285
#define AXIS 286
#define AX 287
#define AY 288
#define AZ 289
#define STRING 290




/* Copy the first part of user declarations.  */
////#line 1 "phantom.y"

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
    //phantomdebug = 1;
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
  
  

/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
//#line 123 "phantom.y"
typedef union YYSTYPE {
  double val;
  double vec[3];
  char str[128];
} YYSTYPE;
/* Line 191 of yacc.c.  */
//#line 283 "phantom.tab.cpp"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 214 of yacc.c.  */
//#line 295 "phantom.tab.cpp"

#if ! defined (yyoverflow) || YYERROR_VERBOSE

# ifndef YYFREE
#  define YYFREE free
# endif
# ifndef YYMALLOC
#  define YYMALLOC malloc
# endif

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   define YYSTACK_ALLOC alloca
#  endif
# else
#  if defined (alloca) || defined (_ALLOCA_H)
#   define YYSTACK_ALLOC alloca
#  else
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
# else
#  if defined (__STDC__) || defined (__cplusplus)
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   define YYSIZE_T size_t
#  endif
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
# endif
#endif /* ! defined (yyoverflow) || YYERROR_VERBOSE */


#if (! defined (yyoverflow) \
     && (! defined (__cplusplus) \
	 || (defined (YYSTYPE_IS_TRIVIAL) && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short) + sizeof (YYSTYPE))				\
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined (__GNUC__) && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  register YYSIZE_T yyi;		\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (0)
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (0)

#endif

#if defined (__STDC__) || defined (__cplusplus)
   typedef signed char yysigned_char;
#else
   typedef short yysigned_char;
#endif

/* YYFINAL -- State number of the termination state. */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   311

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  57
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  47
/* YYNRULES -- Number of rules. */
#define YYNRULES  123
/* YYNRULES -- Number of states. */
#define YYNSTATES  229

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   290

#define YYTRANSLATE(YYX) 						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const unsigned char yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,    43,     2,     2,     2,     2,     2,
      52,    54,    27,    26,    53,    25,     2,    28,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    45,     2,
      55,    24,    56,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    44,     2,    46,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    51,     2,
       2,     2,     2,     2,    50,     2,     2,     2,     2,     2,
      47,    48,    49,    41,     2,    42,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    29,
      30,    31,    32,    33,    34,    35,    36,    37,    38,    39,
      40
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const unsigned short yyprhs[] =
{
       0,     0,     3,     4,     7,    13,    15,    17,    20,    23,
      27,    31,    32,    34,    38,    40,    42,    44,    46,    48,
      50,    52,    54,    56,    58,    60,    62,    64,    66,    68,
      70,    72,    74,    80,    81,    84,    87,    93,    94,    97,
     100,   103,   109,   115,   121,   122,   125,   128,   131,   137,
     138,   141,   144,   147,   150,   156,   157,   160,   162,   164,
     170,   171,   174,   177,   180,   186,   187,   190,   193,   196,
     199,   202,   208,   214,   220,   221,   224,   227,   230,   236,
     237,   240,   243,   246,   249,   255,   261,   267,   268,   271,
     274,   277,   281,   285,   289,   293,   297,   301,   305,   309,
     313,   317,   320,   323,   326,   329,   337,   341,   345,   349,
     353,   357,   361,   366,   371,   373,   378,   383,   388,   392,
     396,   400,   404,   407
};

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const yysigned_char yyrhs[] =
{
      58,     0,    -1,    -1,    58,    59,    -1,    41,    63,    65,
      60,    42,    -1,    61,    -1,    62,    -1,    62,    61,    -1,
      61,    62,    -1,    23,    24,    40,    -1,    19,    24,   103,
      -1,    -1,    64,    -1,    43,    34,    43,    -1,    67,    -1,
      69,    -1,    71,    -1,    72,    -1,    73,    -1,    75,    -1,
      77,    -1,    80,    -1,    82,    -1,    84,    -1,    85,    -1,
      86,    -1,    88,    -1,    90,    -1,    91,    -1,    92,    -1,
      94,    -1,   102,    -1,    44,     3,    45,    68,    46,    -1,
      -1,    68,    95,    -1,    68,    66,    -1,    44,     4,    45,
      70,    46,    -1,    -1,    70,    96,    -1,    70,   100,    -1,
      70,    66,    -1,    44,     5,    45,    74,    46,    -1,    44,
       6,    45,    74,    46,    -1,    44,     7,    45,    74,    46,
      -1,    -1,    74,    95,    -1,    74,    97,    -1,    74,    66,
      -1,    44,     8,    45,    76,    46,    -1,    -1,    76,    95,
      -1,    76,    97,    -1,    76,    99,    -1,    76,    66,    -1,
      44,     9,    45,    78,    46,    -1,    -1,    78,    79,    -1,
      96,    -1,    66,    -1,    44,    10,    45,    81,    46,    -1,
      -1,    81,    96,    -1,    81,   100,    -1,    81,    66,    -1,
      44,    11,    45,    83,    46,    -1,    -1,    83,    99,    -1,
      83,   100,    -1,    83,    96,    -1,    83,    97,    -1,    83,
      66,    -1,    44,    12,    45,    87,    46,    -1,    44,    13,
      45,    87,    46,    -1,    44,    14,    45,    87,    46,    -1,
      -1,    87,    97,    -1,    87,    96,    -1,    87,    66,    -1,
      44,    15,    45,    89,    46,    -1,    -1,    89,    99,    -1,
      89,    97,    -1,    89,    98,    -1,    89,    66,    -1,    44,
      16,    45,    93,    46,    -1,    44,    17,    45,    93,    46,
      -1,    44,    18,    45,    93,    46,    -1,    -1,    89,    97,
      -1,    89,    98,    -1,    89,    66,    -1,    47,    24,   103,
      -1,    48,    24,   103,    -1,    49,    24,   103,    -1,    50,
      24,   103,    -1,    30,    24,   103,    -1,    31,    24,   103,
      -1,    32,    24,   103,    -1,    51,    24,   103,    -1,    33,
      24,   103,    -1,    35,    24,   103,    -1,    36,   101,    -1,
      37,   101,    -1,    38,   101,    -1,    39,   101,    -1,    52,
     103,    53,   103,    53,   103,    54,    -1,    47,    55,   103,
      -1,    47,    56,   103,    -1,    48,    55,   103,    -1,    48,
      56,   103,    -1,    49,    55,   103,    -1,    49,    56,   103,
      -1,    50,   101,    55,   103,    -1,    50,   101,    56,   103,
      -1,    34,    -1,    20,    52,   103,    54,    -1,    21,    52,
     103,    54,    -1,    22,    52,   103,    54,    -1,   103,    26,
     103,    -1,   103,    25,   103,    -1,   103,    27,   103,    -1,
     103,    28,   103,    -1,    25,   103,    -1,    52,   103,    54,
      -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned short yyrline[] =
{
       0,   146,   146,   147,   151,   170,   170,   170,   170,   174,
     178,   181,   182,   186,   189,   190,   191,   192,   193,   194,
     195,   196,   197,   198,   199,   200,   201,   202,   203,   204,
     208,   209,   212,   218,   219,   220,   224,   267,   268,   269,
     270,   274,   282,   290,   297,   298,   299,   300,   304,   315,
     316,   317,   318,   319,   323,   329,   330,   333,   334,   338,
     380,   381,   382,   383,   388,   435,   436,   437,   438,   439,
     440,   444,   452,   460,   466,   467,   468,   469,   473,   476,
     477,   478,   479,   480,   484,   488,   492,   495,   496,   497,
     498,   504,   505,   506,   511,   516,   517,   518,   523,   528,
     529,   533,   537,   538,   539,   543,   552,   559,   566,   573,
     580,   587,   594,   602,   613,   614,   615,   616,   617,   618,
     619,   620,   621,   622
};
#endif

#if YYDEBUG || YYERROR_VERBOSE
/* YYTNME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals. */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "SPHERE", "BOX", "CYLX", "CYLY", "CYLZ",
  "CYL", "ELLIP", "ELLIPFREE", "ELLIPCYL", "ECYLX", "ECYLY", "ECYLZ",
  "CONE", "CONEX", "CONEY", "CONEZ", "RHO", "SIN", "COS", "SQRT",
  "MATERIAL", "'='", "'-'", "'+'", "'*'", "'/'", "NEG", "LENGTHX",
  "LENGTHY", "LENGTHZ", "RADIUS1", "NUMBER", "RADIUS2", "AXIS", "AX", "AY",
  "AZ", "STRING", "'{'", "'}'", "'\"'", "'['", "':'", "']'", "'x'", "'y'",
  "'z'", "'r'", "'l'", "'('", "','", "')'", "'<'", "'>'", "$accept",
  "objectdefs", "objectdef", "dmspec", "materialdef", "densitydef",
  "labeldef", "quotedstring", "volumedef", "commonAssignments",
  "sphereDef", "sphereAssignmentBlock", "boxDef", "boxAssignmentBlock",
  "cylinderXDef", "cylinderYDef", "cylinderZDef",
  "acylinderAssignmentBlock", "cylinderDef", "cylinderAssignmentBlock",
  "ellipsoidDef", "ellipsoidAssignmentBlock",
  "ellipsoidAssignmentStatement", "freeEllipsoidDef",
  "freeEllipsoidAssignmentBlock", "ellipticalCylDef",
  "eCylinderAssignmentBlock", "eCylinderXDef", "eCylinderYDef",
  "eCylinderZDef", "aECylAssignmentBlock", "coneDef",
  "coneAssignmentBlock", "coneXDef", "coneYDef", "coneZDef",
  "aConeAssignmentBlock", "positionAssignment", "radiusAssignment",
  "lengthsAssignment", "lengthAssignment", "radiiAssignment",
  "axisAssignment", "axesAssignment", "vectorExpression", "clipPlanes",
  "expression", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const unsigned short yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,    61,    45,    43,    42,    47,   279,
     280,   281,   282,   283,   284,   285,   286,   287,   288,   289,
     290,   123,   125,    34,    91,    58,    93,   120,   121,   122,
     114,   108,    40,    44,    41,    60,    62
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const unsigned char yyr1[] =
{
       0,    57,    58,    58,    59,    60,    60,    60,    60,    61,
      62,    63,    63,    64,    65,    65,    65,    65,    65,    65,
      65,    65,    65,    65,    65,    65,    65,    65,    65,    65,
      66,    66,    67,    68,    68,    68,    69,    70,    70,    70,
      70,    71,    72,    73,    74,    74,    74,    74,    75,    76,
      76,    76,    76,    76,    77,    78,    78,    79,    79,    80,
      81,    81,    81,    81,    82,    83,    83,    83,    83,    83,
      83,    84,    85,    86,    87,    87,    87,    87,    88,    89,
      89,    89,    89,    89,    90,    91,    92,    93,    93,    93,
      93,    94,    94,    94,    95,    96,    96,    96,    97,    98,
      98,    99,   100,   100,   100,   101,   102,   102,   102,   102,
     102,   102,   102,   102,   103,   103,   103,   103,   103,   103,
     103,   103,   103,   103
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const unsigned char yyr2[] =
{
       0,     2,     0,     2,     5,     1,     1,     2,     2,     3,
       3,     0,     1,     3,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     5,     0,     2,     2,     5,     0,     2,     2,
       2,     5,     5,     5,     0,     2,     2,     2,     5,     0,
       2,     2,     2,     2,     5,     0,     2,     1,     1,     5,
       0,     2,     2,     2,     5,     0,     2,     2,     2,     2,
       2,     5,     5,     5,     0,     2,     2,     2,     5,     0,
       2,     2,     2,     2,     5,     5,     5,     0,     2,     2,
       2,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     2,     2,     2,     2,     7,     3,     3,     3,     3,
       3,     3,     4,     4,     1,     4,     4,     4,     3,     3,
       3,     3,     2,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned char yydefact[] =
{
       2,     0,     1,    11,     3,     0,     0,    12,     0,     0,
       0,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    13,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     5,     6,    33,
      37,    44,    44,    44,    49,    55,    60,    65,    74,    74,
      74,    79,    79,    79,    79,     0,     0,     4,     8,     7,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   114,     0,    10,     9,    32,     0,     0,     0,     0,
      35,    30,    34,    31,     0,     0,     0,     0,     0,     0,
      36,     0,    40,    38,    39,    41,     0,    47,    45,    46,
      42,    43,     0,    48,    53,    50,    51,    52,    54,    58,
      56,    57,    59,    63,    61,    62,    64,    70,    68,    69,
      66,    67,    71,    77,    76,    75,    72,    73,     0,     0,
      78,    83,    81,    82,    80,    83,    81,    82,    84,    85,
      86,     0,     0,     0,   122,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   102,   103,   104,     0,   101,
       0,     0,     0,     0,     0,   123,   119,   118,   120,   121,
      91,   106,   107,    92,   108,   109,    93,   110,   111,    94,
       0,     0,     0,    95,    96,    97,    98,    99,   100,   115,
     116,   117,     0,   112,   113,     0,     0,     0,   105
};

/* YYDEFGOTO[NTERM-NUM]. */
static const short yydefgoto[] =
{
      -1,     1,     4,    46,    47,    48,     6,     7,    10,   117,
      11,    70,    12,    71,    13,    14,    15,    72,    16,    75,
      17,    76,   130,    18,    77,    19,    78,    20,    21,    22,
      79,    23,    83,    24,    25,    26,    84,   101,   118,   144,
     119,   153,   154,   114,   181,   103,    93
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -91
static const short yypact[] =
{
     -91,     1,   -91,   -36,   -91,   -20,    -8,   -91,     7,   244,
      26,   -91,   -91,   -91,   -91,   -91,   -91,   -91,   -91,   -91,
     -91,   -91,   -91,   -91,   -91,   -91,   -91,   -91,   -19,    12,
      16,    21,    54,    68,    75,    86,    92,    93,    97,    98,
      99,   135,   136,   146,    22,    28,    81,   129,   159,   -91,
     -91,   -91,   -91,   -91,   -91,   -91,   -91,   -91,   -91,   -91,
     -91,   -91,   147,   147,   147,   -17,   161,   -91,   -91,   -91,
     259,   137,   218,   247,   253,   227,   168,   158,    78,   103,
     109,   115,   176,   195,   148,   156,   157,   167,   169,   186,
     -17,   -91,   -17,   145,   -91,   -91,   -18,   -15,   -12,   -13,
     -91,   -91,   -91,   -91,   189,   196,   215,   188,   188,   188,
     -91,   188,   -91,   -91,   -91,   -91,   217,   -91,   -91,   -91,
     -91,   -91,   188,   -91,   -91,   -91,   -91,   -91,   -91,   -91,
     -91,   -91,   -91,   -91,   -91,   -91,   -91,   -91,   -91,   -91,
     -91,   -91,   -91,   -91,   -91,   -91,   -91,   -91,   246,   248,
     -91,   -91,   -91,   -91,   -91,   164,   225,   235,   -91,   -91,
     -91,   -17,   -17,   -17,   -91,    -7,   -17,   -17,   -17,   -17,
     -17,   -17,   -17,   -17,   -17,   -17,   -17,   -17,   -17,   -17,
     -17,     3,   -17,   -17,   -17,   -91,   -91,   -91,   -17,   -91,
     -17,   -17,    -3,     2,     6,   -91,    47,    47,   -91,   -91,
     145,   145,   145,   145,   145,   145,   145,   145,   145,   145,
      42,   -17,   -17,   145,   145,   145,   145,   145,   145,   -91,
     -91,   -91,   -17,   145,   145,    77,   -17,    37,   -91
};

/* YYPGOTO[NTERM-NUM].  */
static const short yypgoto[] =
{
     -91,   -91,   -91,   -91,   234,   236,   -91,   -91,   -91,   209,
     -91,   -91,   -91,   -91,   -91,   -91,   -91,    44,   -91,   -91,
     -91,   -91,   -91,   -91,   -91,   -91,   -91,   -91,   -91,   -91,
      52,   -91,   249,   -91,   -91,   -91,    43,   -91,   -60,   -23,
     154,   228,   -62,    41,    70,   -91,   -90
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -91
static const short yytable[] =
{
     164,     2,   165,    87,    88,    89,   170,     5,    90,   173,
     102,   179,   176,   127,     8,   125,   140,    91,   166,   167,
     168,   169,   166,   167,   168,   169,    49,   166,   167,   168,
     169,   166,   167,   168,   169,    92,     9,   171,   172,   180,
     174,   175,     3,   177,   178,    44,    65,   195,   113,    45,
      27,   219,    66,   131,   134,   138,   220,    50,   211,   212,
     221,    51,   166,   167,   168,   169,    52,   166,   167,   168,
     169,   192,   193,   194,   168,   169,   196,   197,   198,   199,
     200,   201,   202,   203,   204,   205,   206,   207,   208,   209,
     210,   228,   213,   214,   215,   222,    73,    74,   216,    53,
     217,   218,   166,   167,   168,   169,    85,    86,   104,   105,
     106,    80,    81,    54,   122,   107,   108,   109,   135,   141,
      55,   223,   224,    67,   136,    96,    97,    98,   111,   116,
     226,    56,   225,   104,   105,   106,   227,    57,    58,   104,
     105,   106,    59,    60,    61,   104,   105,   106,    44,   142,
      96,    97,    98,   111,   116,   146,    96,    97,    98,   111,
     116,   147,    96,    97,    98,   111,   116,   104,   105,   106,
     166,   167,   168,   169,   107,   108,   109,   185,   186,   187,
      62,    63,    45,   110,    96,    97,    98,   111,   104,   105,
     106,    64,   189,   -87,   158,   107,   108,   109,   104,   105,
     106,    94,   159,   160,   132,    96,    97,    98,   111,   148,
     -90,   149,   122,   182,   128,    96,    97,    98,   111,   161,
     183,   162,   150,    96,    97,    98,   111,   116,   148,   126,
     149,   122,   139,   145,   145,   145,   152,   156,   163,   184,
     180,   188,    96,    97,    98,   111,   116,    28,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,    42,    43,   122,   115,    96,    97,    98,    99,   116,
     190,   -88,   191,   123,    96,    97,    98,    99,   116,   100,
     112,   -89,    69,    68,   124,   129,   133,   137,   143,   143,
     143,   151,   155,   120,    96,    97,    98,    99,   116,   121,
      96,    97,    98,    99,   116,    95,    96,    97,    98,    99,
      82,   157
};

static const unsigned char yycheck[] =
{
      90,     0,    92,    20,    21,    22,    24,    43,    25,    24,
      70,    24,    24,    75,    34,    75,    78,    34,    25,    26,
      27,    28,    25,    26,    27,    28,    45,    25,    26,    27,
      28,    25,    26,    27,    28,    52,    44,    55,    56,    52,
      55,    56,    41,    55,    56,    19,    24,    54,    71,    23,
      43,    54,    24,    76,    77,    78,    54,    45,    55,    56,
      54,    45,    25,    26,    27,    28,    45,    25,    26,    27,
      28,   161,   162,   163,    27,    28,   166,   167,   168,   169,
     170,   171,   172,   173,   174,   175,   176,   177,   178,   179,
     180,    54,   182,   183,   184,    53,    52,    53,   188,    45,
     190,   191,    25,    26,    27,    28,    63,    64,    30,    31,
      32,    59,    60,    45,    36,    37,    38,    39,    77,    78,
      45,   211,   212,    42,    46,    47,    48,    49,    50,    51,
      53,    45,   222,    30,    31,    32,   226,    45,    45,    30,
      31,    32,    45,    45,    45,    30,    31,    32,    19,    46,
      47,    48,    49,    50,    51,    46,    47,    48,    49,    50,
      51,    46,    47,    48,    49,    50,    51,    30,    31,    32,
      25,    26,    27,    28,    37,    38,    39,   107,   108,   109,
      45,    45,    23,    46,    47,    48,    49,    50,    30,    31,
      32,    45,   122,    46,    46,    37,    38,    39,    30,    31,
      32,    40,    46,    46,    46,    47,    48,    49,    50,    33,
      46,    35,    36,    24,    46,    47,    48,    49,    50,    52,
      24,    52,    46,    47,    48,    49,    50,    51,    33,    75,
      35,    36,    78,    79,    80,    81,    82,    83,    52,    24,
      52,    24,    47,    48,    49,    50,    51,     3,     4,     5,
       6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    36,    46,    47,    48,    49,    50,    51,
      24,    46,    24,    46,    47,    48,    49,    50,    51,    70,
      71,    46,    48,    47,    75,    76,    77,    78,    79,    80,
      81,    82,    83,    46,    47,    48,    49,    50,    51,    46,
      47,    48,    49,    50,    51,    46,    47,    48,    49,    50,
      61,    83
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
{
       0,    58,     0,    41,    59,    43,    63,    64,    34,    44,
      65,    67,    69,    71,    72,    73,    75,    77,    80,    82,
      84,    85,    86,    88,    90,    91,    92,    43,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    23,    60,    61,    62,    45,
      45,    45,    45,    45,    45,    45,    45,    45,    45,    45,
      45,    45,    45,    45,    45,    24,    24,    42,    62,    61,
      68,    70,    74,    74,    74,    76,    78,    81,    83,    87,
      87,    87,    89,    89,    93,    93,    93,    20,    21,    22,
      25,    34,    52,   103,    40,    46,    47,    48,    49,    50,
      66,    94,    95,   102,    30,    31,    32,    37,    38,    39,
      46,    50,    66,    96,   100,    46,    51,    66,    95,    97,
      46,    46,    36,    46,    66,    95,    97,    99,    46,    66,
      79,    96,    46,    66,    96,   100,    46,    66,    96,    97,
      99,   100,    46,    66,    96,    97,    46,    46,    33,    35,
      46,    66,    97,    98,    99,    66,    97,    98,    46,    46,
      46,    52,    52,    52,   103,   103,    25,    26,    27,    28,
      24,    55,    56,    24,    55,    56,    24,    55,    56,    24,
      52,   101,    24,    24,    24,   101,   101,   101,    24,   101,
      24,    24,   103,   103,   103,    54,   103,   103,   103,   103,
     103,   103,   103,   103,   103,   103,   103,   103,   103,   103,
     103,    55,    56,   103,   103,   103,   103,   103,   103,    54,
      54,    54,    53,   103,   103,   103,    53,   103,    54
};

#if ! defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if ! defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T)
# if defined (__STDC__) || defined (__cplusplus)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# endif
#endif
#if ! defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    { 								\
      yyerror ("syntax error: cannot back up");\
      YYERROR;							\
    }								\
while (0)

#define YYTERROR	1
#define YYERRCODE	256

/* YYLLOC_DEFAULT -- Compute the default location (before the actions
   are run).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)		\
   ((Current).first_line   = (Rhs)[1].first_line,	\
    (Current).first_column = (Rhs)[1].first_column,	\
    (Current).last_line    = (Rhs)[N].last_line,	\
    (Current).last_column  = (Rhs)[N].last_column)
#endif

/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (0)

# define YYDSYMPRINT(Args)			\
do {						\
  if (yydebug)					\
    yysymprint Args;				\
} while (0)

# define YYDSYMPRINTF(Title, Token, Value, Location)		\
do {								\
  if (yydebug)							\
    {								\
      YYFPRINTF (stderr, "%s ", Title);				\
      yysymprint (stderr, 					\
                  Token, Value);	\
      YYFPRINTF (stderr, "\n");					\
    }								\
} while (0)

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_stack_print (short *bottom, short *top)
#else
static void
yy_stack_print (bottom, top)
    short *bottom;
    short *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (/* Nothing. */; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_reduce_print (int yyrule)
#else
static void
yy_reduce_print (yyrule)
    int yyrule;
#endif
{
  int yyi;
  unsigned int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %u), ",
             yyrule - 1, yylno);
  /* Print the symbols being reduced, and their result.  */
  for (yyi = yyprhs[yyrule]; 0 <= yyrhs[yyi]; yyi++)
    YYFPRINTF (stderr, "%s ", yytname [yyrhs[yyi]]);
  YYFPRINTF (stderr, "-> %s\n", yytname [yyr1[yyrule]]);
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (Rule);		\
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YYDSYMPRINT(Args)
# define YYDSYMPRINTF(Title, Token, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   SIZE_MAX < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#if defined (YYMAXDEPTH) && YYMAXDEPTH == 0
# undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined (__GLIBC__) && defined (_STRING_H)
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
#   if defined (__STDC__) || defined (__cplusplus)
yystrlen (const char *yystr)
#   else
yystrlen (yystr)
     const char *yystr;
#   endif
{
  register const char *yys = yystr;

  while (*yys++ != '\0')
    continue;

  return yys - yystr - 1;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined (__GLIBC__) && defined (_STRING_H) && defined (_GNU_SOURCE)
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
#   if defined (__STDC__) || defined (__cplusplus)
yystpcpy (char *yydest, const char *yysrc)
#   else
yystpcpy (yydest, yysrc)
     char *yydest;
     const char *yysrc;
#   endif
{
  register char *yyd = yydest;
  register const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

#endif /* !YYERROR_VERBOSE */



#if YYDEBUG
/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yysymprint (FILE *yyoutput, int yytype, YYSTYPE *yyvaluep)
#else
static void
yysymprint (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  if (yytype < YYNTOKENS)
    {
      YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
# ifdef YYPRINT
      YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
    }
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  switch (yytype)
    {
      default:
        break;
    }
  YYFPRINTF (yyoutput, ")");
}

#endif /* ! YYDEBUG */
/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yydestruct (int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yytype, yyvaluep)
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  switch (yytype)
    {

      default:
        break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM);
# else
int yyparse ();
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM)
# else
int yyparse (YYPARSE_PARAM)
  void *YYPARSE_PARAM;
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  register int yystate;
  register int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  short	yyssa[YYINITDEPTH];
  short *yyss = yyssa;
  register short *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  register YYSTYPE *yyvsp;



#define YYPOPSTACK   (yyvsp--, yyssp--)

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* When reducing, the number of symbols on the RHS of the reduced
     rule.  */
  int yylen;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed. so pushing a state here evens the stacks.
     */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack. Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	short *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow ("parser stack overflow",
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyoverflowlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyoverflowlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	short *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyoverflowlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YYDSYMPRINTF ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */
  YYDPRINTF ((stderr, "Shifting token %s, ", yytname[yytoken]));

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;


  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  yystate = yyn;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 4:
//#line 152 "phantom.y"
    {
  if (top != NULL) {
    //std::cout << "push" << std::endl;
    //std::cout << "x" << x<< std::endl;
    //std::cout << "rho" << rho << std::endl;
    //std::cout << "objectCount" << objectCount << std::endl;
    //std::cout << "ostack size" << ostack.size() << std::endl;
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
;}
    break;

  case 9:
//#line 174 "phantom.y"
#ifndef WIN32
    {material = std::string(strdup(yyvsp[0].str));}
#else
	{material = std::string(_strdup(yyvsp[0].str));}
#endif
    break;

  case 10:
//#line 178 "phantom.y"
    { rho = yyvsp[0].val ;}
    break;

  case 32:
//#line 212 "phantom.y"
    {
  top = new Sphere;
  top->Transform().Scale(Vec3(scale_factor*radius,scale_factor*radius,scale_factor*radius));
;}
    break;

  case 36:
//#line 224 "phantom.y"
    {
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
;}
    break;

  case 41:
//#line 274 "phantom.y"
    {
  top = new Cylinder;
  top->Transform().Scale(Vec3(scale_factor*radius,scale_factor*radius,scale_factor*length));
  top->Transform().RotateY(90);
;}
    break;

  case 42:
//#line 282 "phantom.y"
    {
  top = new Cylinder;
  top->Transform().Scale(Vec3(scale_factor*radius,scale_factor*radius,scale_factor*length));
  top->Transform().RotateX(90);
;}
    break;

  case 43:
//#line 290 "phantom.y"
    {
  top = new Cylinder;
  top->Transform().Scale(Vec3(scale_factor*radius,scale_factor*radius,scale_factor*length));
  top->Transform().RotateZ(90);
;}
    break;

  case 48:
//#line 304 "phantom.y"
    {
  top = new Cylinder;
  top->Transform().Scale(Vec3(scale_factor*radius,scale_factor*radius,scale_factor*length));
   double rotang;
  rotang = -atan2(sqrt(axis_x*axis_x+axis_y*axis_y),axis_z)*180.0/M_PI;
  top->Transform().RotateX(rotang);
  rotang = -atan2(axis_x,axis_y)*180.0/M_PI;
  top->Transform().RotateZ(rotang);
;}
    break;

  case 54:
//#line 323 "phantom.y"
    {
  top = new Sphere;
  top->Transform().Scale(Vec3(scale_factor*lengthx,scale_factor*lengthy,scale_factor*lengthz));
;}
    break;

  case 59:
//#line 338 "phantom.y"
    { 
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
;}
    break;

  case 64:
//#line 388 "phantom.y"
    {
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
;}
    break;

  case 71:
//#line 444 "phantom.y"
    {
  top = new Cylinder();
  top->Transform().Scale(Vec3(scale_factor*lengthz,scale_factor*lengthy,scale_factor*length));
  top->Transform().RotateY(90);
;}
    break;

  case 72:
//#line 452 "phantom.y"
    {
  top = new Cylinder();
  top->Transform().Scale(Vec3(scale_factor*lengthx,scale_factor*lengthz,scale_factor*length));
  top->Transform().RotateX(90);
;}
    break;

  case 73:
//#line 460 "phantom.y"
    {
  top = new Cylinder();
  top->Transform().Scale(Vec3(scale_factor*lengthx,scale_factor*lengthy,scale_factor*length));
;}
    break;

  case 91:
//#line 504 "phantom.y"
    { x = yyvsp[0].val ;}
    break;

  case 92:
//#line 505 "phantom.y"
    { y = yyvsp[0].val ;}
    break;

  case 93:
//#line 506 "phantom.y"
    { z = yyvsp[0].val ;}
    break;

  case 94:
//#line 511 "phantom.y"
    { radius = yyvsp[0].val ;}
    break;

  case 95:
//#line 516 "phantom.y"
    { lengthx = yyvsp[0].val ;}
    break;

  case 96:
//#line 517 "phantom.y"
    { lengthy = yyvsp[0].val ;}
    break;

  case 97:
//#line 518 "phantom.y"
    { lengthz = yyvsp[0].val ;}
    break;

  case 98:
//#line 523 "phantom.y"
    { length = yyvsp[0].val ;}
    break;

  case 99:
//#line 528 "phantom.y"
    { radius1 = yyvsp[0].val ;}
    break;

  case 100:
//#line 529 "phantom.y"
    { radius2 = yyvsp[0].val ;}
    break;

  case 101:
//#line 533 "phantom.y"
    { axis_x = yyvsp[0].vec[0]; axis_y = yyvsp[0].vec[1]; axis_z = yyvsp[0].vec[2]; axis_given = true; ;}
    break;

  case 102:
//#line 537 "phantom.y"
    { ax_x = yyvsp[0].vec[0]; ax_y = yyvsp[0].vec[1]; ax_z = yyvsp[0].vec[2]; ax_given = true; ;}
    break;

  case 103:
//#line 538 "phantom.y"
    { ay_x = yyvsp[0].vec[0]; ay_y = yyvsp[0].vec[1]; ay_z = yyvsp[0].vec[2]; ay_given = true;;}
    break;

  case 104:
//#line 539 "phantom.y"
    { az_x = yyvsp[0].vec[0]; az_y = yyvsp[0].vec[1]; az_z = yyvsp[0].vec[2]; az_given = true;;}
    break;

  case 105:
//#line 544 "phantom.y"
    { yyval.vec[0] = yyvsp[-5].val; yyval.vec[1] = yyvsp[-3].val; yyval.vec[2] = yyvsp[-1].val; ;}
    break;

  case 106:
//#line 552 "phantom.y"
    {
  clipcoefs[clip_count][0] = -1.0;
  clipcoefs[clip_count][1] = 0.0;
  clipcoefs[clip_count][2] = 0.0;
  clipcoefs[clip_count][3] = -scale_factor*yyvsp[0].val;
  clip_count++;
;}
    break;

  case 107:
//#line 559 "phantom.y"
    {
  clipcoefs[clip_count][0] = 1.0;
  clipcoefs[clip_count][1] = 0.0;
  clipcoefs[clip_count][2] = 0.0;
  clipcoefs[clip_count][3] = scale_factor*yyvsp[0].val;
  clip_count++;
;}
    break;

  case 108:
//#line 566 "phantom.y"
    {
  clipcoefs[clip_count][0] = 0.0;
  clipcoefs[clip_count][1] = -1.0;
  clipcoefs[clip_count][2] = 0.0;
  clipcoefs[clip_count][3] = -scale_factor*yyvsp[0].val;
  clip_count++;
;}
    break;

  case 109:
//#line 573 "phantom.y"
    {
  clipcoefs[clip_count][0] = 0.0;
  clipcoefs[clip_count][1] = 1.0;
  clipcoefs[clip_count][2] = 0.0;
  clipcoefs[clip_count][3] = scale_factor*yyvsp[0].val;
  clip_count++;
;}
    break;

  case 110:
//#line 580 "phantom.y"
    {
  clipcoefs[clip_count][0] = 0.0;
  clipcoefs[clip_count][1] = 0.0;
  clipcoefs[clip_count][2] = -1.0;
  clipcoefs[clip_count][3] = -scale_factor*yyvsp[0].val;
  clip_count++;
;}
    break;

  case 111:
//#line 587 "phantom.y"
    {
  clipcoefs[clip_count][0] = 0.0;
  clipcoefs[clip_count][1] = 0.0;
  clipcoefs[clip_count][2] = 1.0;
  clipcoefs[clip_count][3] = scale_factor*yyvsp[0].val;
  clip_count++;
;}
    break;

  case 112:
//#line 594 "phantom.y"
    {
  double norm = sqrt(yyvsp[-2].vec[0]*yyvsp[-2].vec[0] + yyvsp[-2].vec[1]*yyvsp[-2].vec[1] + yyvsp[-2].vec[2]*yyvsp[-2].vec[2]);
  clipcoefs[clip_count][0] = -yyvsp[-2].vec[0]/norm;
  clipcoefs[clip_count][1] = -yyvsp[-2].vec[1]/norm;
  clipcoefs[clip_count][2] = -yyvsp[-2].vec[2]/norm;
  clipcoefs[clip_count][3] = -scale_factor*yyvsp[0].val;
  clip_count++;
;}
    break;

  case 113:
//#line 602 "phantom.y"
    {
  double norm = sqrt(yyvsp[-2].vec[0]*yyvsp[-2].vec[0] + yyvsp[-2].vec[1]*yyvsp[-2].vec[1] + yyvsp[-2].vec[2]*yyvsp[-2].vec[2]);
  clipcoefs[clip_count][0] = yyvsp[-2].vec[0]/norm;
  clipcoefs[clip_count][1] = yyvsp[-2].vec[1]/norm;
  clipcoefs[clip_count][2] = yyvsp[-2].vec[2]/norm;
  clipcoefs[clip_count][3] = scale_factor*yyvsp[0].val;
  clip_count++;
;}
    break;

  case 115:
//#line 614 "phantom.y"
    { yyval.val = sin(yyvsp[-1].val*M_PI/180.0) ;}
    break;

  case 116:
//#line 615 "phantom.y"
    { yyval.val = cos(yyvsp[-1].val*M_PI/180.0) ;}
    break;

  case 117:
//#line 616 "phantom.y"
    { yyval.val = sqrt(yyvsp[-1].val) ;}
    break;

  case 118:
//#line 617 "phantom.y"
    { yyval.val = yyvsp[-2].val + yyvsp[0].val ;}
    break;

  case 119:
//#line 618 "phantom.y"
    { yyval.val = yyvsp[-2].val - yyvsp[0].val ;}
    break;

  case 120:
//#line 619 "phantom.y"
    { yyval.val = yyvsp[-2].val * yyvsp[0].val ;}
    break;

  case 121:
//#line 620 "phantom.y"
    { yyval.val = yyvsp[-2].val / yyvsp[0].val ;}
    break;

  case 122:
//#line 621 "phantom.y"
    { yyval.val = -yyvsp[0].val ;}
    break;

  case 123:
//#line 622 "phantom.y"
    { yyval.val = yyvsp[-1].val ;}
    break;


    }

/* Line 1000 of yacc.c.  */
//#line 1877 "phantom.tab.cpp"

  yyvsp -= yylen;
  yyssp -= yylen;


  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;
  

/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (YYPACT_NINF < yyn && yyn < YYLAST)
	{
	  YYSIZE_T yysize = 0;
	  int yytype = YYTRANSLATE (yychar);
	  const char* yyprefix;
	  char *yymsg;
	  int yyx;

	  /* Start YYX at -YYN if negative to avoid negative indexes in
	     YYCHECK.  */
	  int yyxbegin = yyn < 0 ? -yyn : 0;

	  /* Stay within bounds of both yycheck and yytname.  */
	  int yychecklim = YYLAST - yyn;
	  int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
	  int yycount = 0;

	  yyprefix = ", expecting ";
	  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	      {
		yysize += yystrlen (yyprefix) + yystrlen (yytname [yyx]);
		yycount += 1;
		if (yycount == 5)
		  {
		    yysize = 0;
		    break;
		  }
	      }
	  yysize += (sizeof ("syntax error, unexpected ")
		     + yystrlen (yytname[yytype]));
	  yymsg = (char *) YYSTACK_ALLOC (yysize);
	  if (yymsg != 0)
	    {
	      char *yyp = yystpcpy (yymsg, "syntax error, unexpected ");
	      yyp = yystpcpy (yyp, yytname[yytype]);

	      if (yycount < 5)
		{
		  yyprefix = ", expecting ";
		  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
		    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
		      {
			yyp = yystpcpy (yyp, yyprefix);
			yyp = yystpcpy (yyp, yytname[yyx]);
			yyprefix = " or ";
		      }
		}
	      yyerror (yymsg);
	      YYSTACK_FREE (yymsg);
	    }
	  else
	    yyerror ("syntax error; also virtual memory exhausted");
	}
      else
#endif /* YYERROR_VERBOSE */
	yyerror ("syntax error");
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* If at end of input, pop the error token,
	     then the rest of the stack, then return failure.  */
	  if (yychar == YYEOF)
	     for (;;)
	       {
		 YYPOPSTACK;
		 if (yyssp == yyss)
		   YYABORT;
		 YYDSYMPRINTF ("Error: popping", yystos[*yyssp], yyvsp, yylsp);
		 yydestruct (yystos[*yyssp], yyvsp);
	       }
        }
      else
	{
	  YYDSYMPRINTF ("Error: discarding", yytoken, &yylval, &yylloc);
	  yydestruct (yytoken, &yylval);
	  yychar = YYEMPTY;

	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

#ifdef __GNUC__
  /* Pacify GCC when the user code never invokes YYERROR and the label
     yyerrorlab therefore never appears in user code.  */
  if (0)
     goto yyerrorlab;
#endif

  yyvsp -= yylen;
  yyssp -= yylen;
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;

      YYDSYMPRINTF ("Error: popping", yystos[*yyssp], yyvsp, yylsp);
      yydestruct (yystos[yystate], yyvsp);
      YYPOPSTACK;
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  YYDPRINTF ((stderr, "Shifting error token, "));

  *++yyvsp = yylval;


  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*----------------------------------------------.
| yyoverflowlab -- parser overflow comes here.  |
`----------------------------------------------*/
yyoverflowlab:
  yyerror ("parser stack overflow");
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}





