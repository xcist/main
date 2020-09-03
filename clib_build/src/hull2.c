// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

/*

  hull2 computes the convex hull of a finite set of points in a plane

  It takes the 2D coordinates as inputs and produces a list of index values of the points on the boundary of the hull in counter clockwise order.

 */

struct pnt
{
  double x;
  double y;
  double tan_angle;
  short index;
};

int compare_pts (const void *a, const void *b)

{
  const struct pnt *da = (const struct pnt *) a;
  const struct pnt *db = (const struct pnt *) b;
  
  return (da[0].tan_angle > db[0].tan_angle) - (da[0].tan_angle < db[0].tan_angle);
}

int hull2(double* x,double* y,short pts,short* list)

{
  short min_index=0,i,list_length=2;
  double min_y=y[0];
  struct pnt *points = NULL;
  const short one_based = 1;

  for(i=1;i<pts;i++)
    if (y[i]<min_y) 
      {min_y=y[i];min_index=i;}

  points = malloc(sizeof(struct pnt)*pts);
  for(i=0;i<pts;i++)
    {
      points[i].x = x[i];
      points[i].y = y[i];
      points[i].index = i;
      points[i].tan_angle = -(x[i] - x[min_index]) / (y[i] - y[min_index]);
    }
  points[min_index].tan_angle = -1e300;

  //for(i=0;i<pts;i++) printf("angle: %1.7lf  index: %d\n\r",points[i].tan_angle,points[i].index);

  qsort(points,pts,sizeof(struct pnt),compare_pts);

  //for(i=0;i<pts;i++) printf("angle: %1.7lf  index: %d\n\r",points[i].tan_angle,points[i].index);

  list[0]=0;list[1]=1;
  
  for(i=2;i<pts;i++)
    {
      while ((list_length >= 2) && (((points[list[list_length-1]].x - points[list[list_length-2]].x) * (points[i].y - points[list[list_length-2]].y) - (points[list[list_length-1]].y - points[list[list_length-2]].y) * (points[i].x - points[list[list_length-2]].x)) <= 0))
	  list_length--;
      list[list_length] = i;
      list_length++;
    }
  for(i=0;i<list_length;i++)
    list[i] = points[list[i]].index + one_based;

  free(points);
  return (int) list_length;

}


// gcc -shared -fPIC -O3 -mfpmath=sse -msse -o hull2.so hull2.c
// gcc -shared -fPIC -O3 -mfpmath=sse -msse -m32 -o hull2.so hull2.c


