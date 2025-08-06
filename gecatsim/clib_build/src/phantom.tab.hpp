/* -----------------------------------------------------------------------
*   Program Name: phantom.tab.hpp        
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




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 123 "phantom.y"
typedef union YYSTYPE {
  double val;
  double vec[3];
  char str[128];
} YYSTYPE;
/* Line 1275 of yacc.c.  */
#line 113 "phantom.tab.hpp"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE phantomlval;




