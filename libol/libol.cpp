/*
        OpenLase - a realtime laser graphics toolkit

Copyright (C) 2009-2011 Hector Martin "marcan" <hector@marcansoft.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 2.1 or version 3.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <QtDebug>

#include "libol.h"
//#include <jack/jack.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
//#include <unistd.h>


#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"


/////////////////////////////////////////////////////////////////////////
//                         Data Structures
/////////////////////////////////////////////////////////////////////////

/*
    struct Object, representing an primitive Object in openlase, similar to openGL primitives
    used in struct Frame and struct DrawState

    int pointcnt, point count
    Point * points, a series of points for a given object
*/

typedef struct {
	int pointcnt;
	Point *points;
} Object;

/*
    struct Frame, current working frame
    used in render engine as working frame buffer

    int objcnt, current object count in working frame
    int objmax, maxium object allowable for allocation
    Object *objects, pointer to the series of objects
    int psmax, points max, maximum number of points allowed for this working frame
    int psnext, offset to the last+1 element in the array of points
    Point *point, pointer to array of points

*/

typedef struct {
	int objcnt;
	int objmax;
	Object *objects;
	int psmax;
	int psnext;
	Point *points;
} Frame;


/*
    Pointer to array of RenderFrames
    used as buffers so frames can be pre-rendered
*/

static RenderedFrame *frames;


/*
    Specifies last frame's info pertaining to last rendering operation
*/

static OLFrameInfo last_info;

/*
    Allocation of working frame
*/

static Frame wframe;

/*
    Working data structure used in drawing operations. Keeps track of current object, last point, last slope, points of concern, primitive to draw, state and number of points
    Essential tool for drawing an object
*/
typedef struct {
	Object *curobj;
	Point last_point;
	Point last_slope;
	Point c1, c2;
	int prim;
	int state;
	int points;
} DrawState;

/*
    Allocation of DrawState
*/
static DrawState dstate;

/*
    Params that specify drawing operations
*/
static OLRenderParams params;

static Point last_render_point;

/*
   hardcoded values
*/

static volatile int cwbuf = 0;
static int fbufs = 1;

/*
    Bounding Box for points {[-1,1], [-1, 1]}
*/
float bbox[2][2];

/*
    Defines depth of 2D transformation matrix stack
*/
#define MTX_STACK_DEPTH 16

/*
    Index for 2D transformation matrix stack
*/
int mtx2dp = 0;

/*
    Allocation for 2D transformation matrix stack
*/
float mtx2ds[MTX_STACK_DEPTH][3][3];

/*
    Current working 2D transformation matrix
*/
float mtx2d[3][3];

/*
    Index for the 3D transfomration matrix stack
*/
int mtx3dp = 0;

/*
    Allocation of 3D transformation matrix stack
*/
float mtx3ds[MTX_STACK_DEPTH][4][4];

/*
    Current working 3D transformation matrix
*/
float mtx3d[4][4];

/*
    Index for the color stack
*/
int coldp = 0;

/*
    Allocation for the color stack
*/
uint32_t cols[MTX_STACK_DEPTH];


/*
    Current color
*/
uint32_t curcol;

/*
    Macro for making a struct point
*/
//#define POINT(x, y, color) ((Point){x,y,color})
#define POINT(x,y, color) Point{x, y, color};

//static inline Point make_point_struct(float x, float y, uint32_t color){
//    Point p = {x, y, color};
//    return p;
//}

/*
    Various shaders
*/
ShaderFunc vpreshader;
ShaderFunc vshader;
Shader3Func v3shader;
ShaderFunc pshader;

/*
    Audio Callback function
*/
AudioCallbackFunc audiocb;

/*
    Logging Callback function
*/
LogCallbackFunc log_cb;


/*****************************************************
 *
 * Start of internal operations
 *
 * **************************************************/

/*
    Color multiplication
*/
static uint32_t colmul(uint32_t a, uint32_t b)
{
    uint32_t out = 0;
    out |= ((a&0xff)*(b&0xff)) / 255;
    out |= ((((a>>8)&0xff)*((b>>8)&0xff)) / 255)<<8;
    out |= ((((a>>16)&0xff)*((b>>16)&0xff)) / 255)<<16;
    return out;
}

/*
    Returns pointer to first available location in points array in the working frame.
    Pushes back psnext index by count units so count points are reserved for next operation.
*/
static Point *ps_alloc(int count)
{
	Point *ret;
	if ((count + wframe.psnext) > wframe.psmax) {
        olLog("Point buffer overflow (temp): need %d points, have %d\n", count + wframe.psnext, wframe.psmax);
        return 0;
//		exit(1);
	}
	ret = wframe.points + wframe.psnext;
	wframe.psnext += count;
	return ret;
}



/*****************************************************
 *
 * End  of operations
 *
 * ****************************************************/

/////////////////////////////////////////////////////////////////////////
//                Libol Module Operations
/////////////////////////////////////////////////////////////////////////



/*
     initialization function
*/
int olInit(int max_points)
{
	int i;

    /*
        zero out dstate and last_render_point
    */
    memset(&dstate, 0, sizeof(dstate));
	memset(&last_render_point, 0, sizeof(last_render_point));



    /*
        initialize working frame.
        initial object max is 16
        allocate array of 16 objects
        max points is taken from olInit argument
        then allocate array of psmax points
    */
    memset(&wframe, 0, sizeof(Frame));
	wframe.objmax = 16;
    wframe.objects = (Object *) malloc(wframe.objmax * sizeof(Object));
	wframe.psmax = max_points;
    wframe.points = (Point *) malloc(wframe.psmax * sizeof(Point));


    /*
        allocate the rendered frame buffer
    */
    frames = (RenderedFrame*) malloc(fbufs * sizeof(RenderedFrame));
    for (i=0; i<fbufs; i++) {
		memset(&frames[i], 0, sizeof(RenderedFrame));
		frames[i].pmax = max_points;
        frames[i].points = (Point *) malloc(frames[i].pmax * sizeof(Point));
        frames[i].audio_l = (float*) malloc(frames[i].pmax * sizeof(float));
        frames[i].audio_r = (float*) malloc(frames[i].pmax * sizeof(float));
    }

    /*
        This resets 2d transformation matrix to identity matrix
    */
    olLoadIdentity();
	for(i=0; i<MTX_STACK_DEPTH; i++)
		olPushMatrix();
	mtx2dp = 0;

    /*
        This resets 3d transformation matrix to identity matrix
     */
    olLoadIdentity3();
	for(i=0; i<MTX_STACK_DEPTH; i++)
		olPushMatrix3();
	mtx3dp = 0;

    /*
        Resets color matrix
    */
    olResetColor();
	for(i=0; i<MTX_STACK_DEPTH; i++)
		olPushColor();
	coldp = 0;

    /*
        Initialize bounding box
    */
    bbox[0][0] = -1;
	bbox[0][1] = -1;
	bbox[1][0] = 1;
	bbox[1][1] = 1;

    /*
        Resets shades and audio call back
    */
    vpreshader = NULL;
	vshader = NULL;
	v3shader = NULL;
	pshader = NULL;
	audiocb = NULL;

	return 0;
}

void olSetRenderParams(OLRenderParams *sp)
{
	params = *sp;
}

void olGetRenderParams(OLRenderParams *sp)
{
	*sp = params;
}


/*
    Starts drawing an object
*/
void olBegin(int prim)
{
    /* If there is currently a drawing object, return and do nothing. dstate.curobj needs to be null in order to start to draw a new object. */
    if (dstate.curobj)
		return;
    /* If working frame's object count exceeds maximum allowed, allocate twice as many*/
    if (wframe.objmax == wframe.objcnt) {
		wframe.objmax *= 2;
        wframe.objects = (Object*) realloc(wframe.objects, wframe.objmax * sizeof(Object));
	}
    /* Set the dstate.curobj to the next object on working frame and zero it */
    dstate.curobj = wframe.objects + wframe.objcnt;
	memset(dstate.curobj, 0, sizeof(Object));
    /* Set the points pointer of current obj to the next points pointer in working frame*/
    dstate.curobj->points = wframe.points + wframe.psnext;
    /* Specifies which primtive to draw */
    dstate.prim = prim;
    /* @question: what is state used for? */
    dstate.state = 0;
    /* @question: what is points used for?*/
    dstate.points = 0;
}


/////////////////////////////////////////////////////////////////////////
//                         Graphic Operations I
/////////////////////////////////////////////////////////////////////////

/*
    checks if a point is near another.
    returns 1 if distance between two points are <= params.snap. 0 otherwise.
*/
static int near(Point a, Point b)
{
	float dx = a.x - b.x;
	float dy = a.y - b.y;
	return sqrtf(dx*dx+dy*dy) <= params.snap;
}

/*
    Adds one point to current object, incrementing point count.
*/
static void addpoint(float x, float y, uint32_t color)
{
    if (dstate.curobj->pointcnt >= 30000) {
        qDebug() << "object has many points" << dstate.curobj->pointcnt;
        return;
    }
    Point *pnt = ps_alloc(1);
	pnt->x = x;
	pnt->y = y;
	pnt->color = color;
	dstate.curobj->pointcnt++;
}

/// Get dwell points number for system current state
static int get_dwell(float x, float y)
{
    // first object point was drawn?
    if (dstate.points == 1) {
        // return settings start dwell
        return params.start_dwell;
	} else {
        // gonna calculate andgle between previous and current drawn segments
        //////////////////////////////////////////////////////////////////////

        // get last drawn point
        float ex = dstate.last_point.x;
		float ey = dstate.last_point.y;
        // get drawn point before last (previous drawn segment)
        float ecx = dstate.last_slope.x;
		float ecy = dstate.last_slope.y;
        // get last drawn point
        float sx = ex;
		float sy = ey;
        // get next drawn point (current drawn segment)
        float scx = x;
		float scy = y;
        // coordinates of previous drawn segment
        float dex = ecx-ex;
		float dey = ecy-ey;
        // coordinates of current drawn segment
        float dsx = sx-scx;
		float dsy = sy-scy;
        // get dot product of two segments
        float dot = dex*dsx + dey*dsy;
        // get product of drawn segments lenghts
        float lens = sqrtf(dex*dex+dey*dey) * sqrtf(dsx*dsx+dsy*dsy);
		//olLog("%f,%f -> %f,%f -> %f,%f\n", ecx,ecy,ex,ey,x,y);
        // TODO: was "if (lens < 0.0001f) {" before
        // check if one of the segments has zero lenght
        if (lens == 0) {
            // use corner dwell value
            //olLog("deg cor\n");
			return params.corner_dwell;
		} else {
            // get angle cosine
            dot = dot / lens;
            // check if angle cosine more than settings value (angle is less!!!)
            if (dot > params.curve_angle) {
                // use curve dwell value
                //olLog("curve\n");
				return params.curve_dwell;
			} else {
                // use corner dwell value
                //olLog("cor\n");
				return params.corner_dwell;
			}
		}
	}
}

/// Draw line (internal)
static void line_to(float x, float y, uint32_t color)
{
    int dwell, i;
	//olLog("points: %d %d\n", dstate.points, dstate.curobj->pointcnt	);
    // if it's first object point
	if (dstate.points == 0) {
        // add point
        addpoint(x,y,color);
        // increase points counter
        dstate.points++;
        // save point as last drawn point
        dstate.last_point = POINT(x,y,color);
		return;
	}
    // get dwell points number
    dwell = get_dwell(x, y);
    // drawn last point dwell value times
    Point last = dstate.last_point;
	for (i=0; i<dwell; i++)
		addpoint(last.x,last.y,last.color);
    // gonna calculate how many points we need to draw to the target point according to settings speed
    /////////////////////////////////////////////////////////////////////////////////////////////

    // get X-distance to target point
    float dx = x - last.x;
    // get Y-distance to target point
    float dy = y - last.y;
    // get maximum distance
    float distance = fmaxf(fabsf(dx),fabsf(dy));
    // get points number
    int points = ceilf(distance/params.on_speed);
    // draw intermediate points
    for (i=1; i<=points; i++) {
		addpoint(last.x + (dx/(float)points) * i, last.y + (dy/(float)points) * i, color);
	}
    //if (points > 1) qDebug() << "points is " << points;
    // save last point as slope point
    dstate.last_slope = dstate.last_point;
    // save target point as last point
    dstate.last_point = POINT(x,y,color);
    // increase points counter
    dstate.points++;
}

// Recursive Bezier interpolation
static void recurse_bezier(float x1, float y1, float x2, float y2, float x3, float y3, uint32_t color, int depth)
{
    // get last drawn point
    float x0 = dstate.last_point.x;
	float y0 = dstate.last_point.y;
    // need subdivide flag
    int subdivide = 0;

    // check recursive depth
    // TODO: was 50 before
	if (depth > 100) {
        // interrupt recursion!
        olLog("Bezier recurse error: %f,%f %f,%f %f,%f %f,%f\n", x0, y0, x1, y1, x2, y2, x3, y3);
		return;
	}
	
    // get X-distance between last and target points
    float dx = x3-x0;
    // get Y-distance between last and target points
    float dy = y3-y0;
    // get maximum distance
    float distance = fmaxf(fabsf(dx),fabsf(dy));
    // compare distance is settings speed
    if (distance > params.on_speed) {
        // if distance is more - divide
        subdivide = 1;
	} else {
        // if distance is less - check flatness
        ////////////////////////////////////////

        // get second point dispersion
        float ux = (3.0*x1 - 2.0*x0 - x3); ux = ux * ux;
		float uy = (3.0*y1 - 2.0*y0 - y3); uy = uy * uy;
        // get third point dispersion
        float vx = (3.0*x2 - 2.0*x3 - x0); vx = vx * vx;
		float vy = (3.0*y2 - 2.0*y3 - y0); vy = vy * vy;
        // get maximum X-dispersion
        if (ux < vx)
			ux = vx;
        // get maximum Y-dispersion
        if (uy < vy)
			uy = vy;
        // compare dispersions sum to settings flatness
        if ((ux+uy) > params.flatness)
			subdivide = 1;
	}

    // TODO: was added by Wicked Lasers??
    if (depth > 5) subdivide = false; // shortcut

    // check we need to divide current spline in two subsplines
    if (subdivide) {
        // interpolate spline with parameter 0.5
        //
        //                                     (m!)
        //                           (a1) (a2) (mc) (b1) (b2)
        // x0___x1___x2___x3 ==> x0________x1________x2________x3
        //////////////////////////////////////////////////////////

        // generate medium point for second and third points
        //de Casteljau at t=0.5
		float mcx = (x1 + x2) * 0.5;
		float mcy = (y1 + y2) * 0.5;
        // generate medium point for first and second points
        float ax1 = (x0 + x1) * 0.5;
		float ay1 = (y0 + y1) * 0.5;
        // generate interpolated point instead x1
        float ax2 = (ax1 + mcx) * 0.5;
		float ay2 = (ay1 + mcy) * 0.5;
        // generated medium point for third and fourth points
        float bx2 = (x2 + x3) * 0.5;
		float by2 = (y2 + y3) * 0.5;
        // generate interpolated point instead x2
        float bx1 = (bx2 + mcx) * 0.5;
		float by1 = (by2 + mcy) * 0.5;
        // generate interpolated point instead mc
        float xm = (ax2 + bx1) * 0.5;
		float ym = (ay2 + by1) * 0.5;
        // interpolate recursively the first part - a1___a2___m
        recurse_bezier(ax1, ay1, ax2, ay2, xm, ym, color, depth+1);
        // interpolate recursively the second part - b1___b2___x3
        recurse_bezier(bx1, by1, bx2, by2, x3, y3, color, depth+1);
	} else {
        // add target point
        addpoint(x3, y3, color);
        // set target point as last drawn point
        dstate.last_point = POINT(x3,y3,color);
	}
}

// Draw bezier spline
static void bezier_to(float x, float y, uint32_t color)
{
    // if it's first object point
    int dwell, i;

	if (dstate.points == 0) {
        // add point
        addpoint(x,y,color);
        // increase points counter
        dstate.points++;
        // set point as last drawn point
        dstate.last_point = POINT(x,y,color);
		return;
	}

    // we need at lest three points (c1, c2 and target point) for interpolation!!!
    ///////////////////////////////////////////////////////////////////////////////

    // check current state
    switch(dstate.state) {
		case 0:
            // we draw first point - save it to c1
            dstate.c1 = POINT(x,y,color);
            // change current state
            dstate.state++;
			return;
		case 1:
            // we draw second point - save it to c2
            dstate.c2 = POINT(x,y,color);
            // change current state
            dstate.state++;
			return;
		case 2:
            // we draw third or more point - interpolation!
            break;
	}

    // check if c1 point is near to last drawn point
    if (near(dstate.last_point, dstate.c1))
        // use c2 point to calculate dwell value
        dwell = get_dwell(dstate.c2.x, dstate.c2.y);
	else
        // use c1 point to calculate dwell value
        dwell = get_dwell(dstate.c1.x, dstate.c1.y);

    // draw last point dwell number times
    Point last = dstate.last_point;
	for (i=0; i<dwell; i++)
		addpoint(last.x,last.y,last.color);

    // start interpolation recursion
    recurse_bezier(dstate.c1.x, dstate.c1.y, dstate.c2.x, dstate.c2.y, x, y, color, 0);

    // save target point as last drawn point
    dstate.last_point = POINT(x,y,color);
    // check if target point is near to c2 point
    if (near(dstate.c2, dstate.last_point))
        // set c1 point as slope point
        dstate.last_slope = dstate.c1;
	else
        // set c2 point as slope point
        dstate.last_slope = dstate.c2;
    // increase points counter
    dstate.points++;
    // reset current state
    dstate.state = 0;
}


/// Draw point (internal)
static void point_to(float x, float y, uint32_t color)
{
    int i;
    // add target point
	addpoint(x,y,color);
    // check if target point is first drawn object point
    if (dstate.points == 0)
        // repeate target point start_dwell value times
        for (i=0; i<params.start_dwell; i++)
			addpoint(x,y,color);

    // increase drawn points counter
    dstate.points++;
	return;
}

// Draw vertex according to current drawing mode
void olVertex(float x, float y, uint32_t color)
{
    // check args
    if (x != x || y != y)
    {
        qDebug() << " x " << x << " y " << y;
        return;
    }

    // check if current object is presented (olBegin() was called)
    if (!dstate.curobj)
		return;

	float nx, ny;

    // check if vertex pre-shader is set
    if(vpreshader)
        // update point
        vpreshader(&x, &y, &color);

    // get blended color (between current and target)
    color = colmul(color,curcol);

    // apply transform 2D matrix
    nx = mtx2d[0][0] * x + mtx2d[0][1] * y + mtx2d[0][2];
	ny = mtx2d[1][0] * x + mtx2d[1][1] * y + mtx2d[1][2];

    // check valid
    if (nx != nx || ny != ny)
    {
        qDebug() << " nx " << nx << " ny " << ny;
        return;
    }

    // check if vertex shader is set
    if(vshader)
        // update point
        vshader(&nx, &ny, &color);

    // check current drawing mode
    switch (dstate.prim) {
		case OL_LINESTRIP:
            // draw line
            line_to(nx,ny,color);
			break;
		case OL_BEZIERSTRIP:
            // draw bezier spline
            bezier_to(nx,ny,color);
			break;
		case OL_POINTS:
            // draw point
            point_to(nx,ny,color);
			break;
	}
}

/// Stop drawing
void olEnd(void)
{
    // check current object
    int i;
	if (!dstate.curobj)
		return;
    // check current object has at least two points
    if (dstate.points < 2) {
        //qDebug() << "nff jfo obj";
        dstate.curobj = NULL;
		return;
	}
    //qDebug() << "end with pts" << dstate.curobj->pointcnt;
    // get current object last point
    Point *last = dstate.curobj->points + dstate.curobj->pointcnt - 1;
    // repeate last point end_dwell times
    for (i=0; i<params.end_dwell; i++)
		addpoint(last->x,last->y,last->color);

    // check pixel shader
    if(pshader) {
        // apply pixel shader to every current object point
        for (i=0; i<dstate.curobj->pointcnt; i++) {
			pshader(&dstate.curobj->points[i].x, &dstate.curobj->points[i].y, &dstate.curobj->points[i].color);
		}
	}

    // check if at leaast one point is in [-1 : +1] for X-axis and for Y-axis
    int nl=0,nr=0,nu=0,nd=0;
	for (i=0; i<dstate.curobj->pointcnt; i++) {
        // skip black point
        if (!dstate.curobj->points[i].color)
			continue;
        // check left bound
        if (dstate.curobj->points[i].x > -1)
			nl = 1;
        // check right bound
        if (dstate.curobj->points[i].x < 1)
			nr = 1;
        // check down bound
        if (dstate.curobj->points[i].y > -1)
			nd = 1;
        // check up bound
        if (dstate.curobj->points[i].y < 1)
			nu = 1;

        // break if point was found
        if (nl && nr && nu && nd)
			break;
	}
    if (nl && nr && nu && nd)
    // increase visible objects counter
		wframe.objcnt++;
    // reset current object pointer
    dstate.curobj = NULL;
}

/// Check necessary space for current object points buffer
static void chkpts(int count)
{
	if (frames[cwbuf].pnext + count > frames[cwbuf].pmax) {
        qDebug() << "checkpts - point overflow";
        olLog("Point buffer overflow (final): need %d points, have %d\n",
				count + frames[cwbuf].pnext, frames[cwbuf].pmax);
        // VERY BAD CODE!!! NOT APPLICABLE FOR APPLE APPSTORE!!!
//		exit(1);
	}
}

/// Add render point for the render buffer
static void addrndpoint(float x, float y, uint32_t color)
{
    if (frames[0].pnext >= frames[0].pmax) { qDebug() << "point overflow"; return; }
    frames[cwbuf].points[frames[cwbuf].pnext].x = x;
	frames[cwbuf].points[frames[cwbuf].pnext].y = y;
	frames[cwbuf].points[frames[cwbuf].pnext].color = color;
	frames[cwbuf].pnext++;
}

/// Convert object points to render points
static void render_object(Object *obj)
{
    // we need to check all object points and skip points are not in the bounding box

    int i,j;
    // first object point pointer
	Point *start = &obj->points[0];
    // last object point pointer
    Point *end = &obj->points[obj->pointcnt-1];
    // get distance from the last drawn point to the first object point
    float dx = start->x - last_render_point.x;
	float dy = start->y - last_render_point.y;
	float distance = fmaxf(fabsf(dx),fabsf(dy));
    // calculate need points to track laser
    int points = ceilf(distance/params.off_speed);
    // check if render buffer has necessary size
    chkpts(2 * (obj->pointcnt + params.start_wait + params.end_wait + points));

    /*qDebug() << params.start_dwell << " "
        << params.start_wait << " "
        << params.end_wait << " "
        << params.end_dwell << " ";
        */

    // last point outside
    Point *out_start = NULL;
    // skip end wait pointers flag
    int skip_out_start_wait = 0;

	Point *ip = obj->points;
    // check if at least one point is in bounding box and has no block color
    for (i=0; i<obj->pointcnt; i++, ip++) {
		if (ip->x < bbox[0][0] || ip->x > bbox[1][0] ||
			ip->y < bbox[0][1] || ip->y > bbox[1][1])
			continue;
		if (ip->color != C_BLACK)
			break;
	}
    // return if points was not found
    if (i == obj->pointcnt) // null object
		return;

    // check if first point is out of bounding box
    if (start->x < bbox[0][0] || start->x > bbox[1][0] ||
		start->y < bbox[0][1] || start->y > bbox[1][1]) {
        // set last drawn point as last point outside
        out_start = &last_render_point;
        // set skip end wait points flag
        skip_out_start_wait = 1;
	} else if (distance > params.snap) {
        // add tracking render points
        for (i=0; i<points; i++) {
			addrndpoint(last_render_point.x + (dx/(float)points) * i,
						last_render_point.y + (dy/(float)points) * i,
						C_BLACK);
		}
        // add start wait render points
        for (i=0; i<params.start_wait; i++) {
			addrndpoint(start->x, start->y, C_BLACK);
		}
	}
    // pointer to the next render point
    if (frames[0].pnext >= frames[0].pmax) { qDebug() << "point overflow"; return; }// else { qDebug() << "asdf \n" << frames[0].pnext; }
    Point *op = &frames[cwbuf].points[frames[cwbuf].pnext];
    // pointer to the object points
    ip = obj->points;
    // check all object points
    for (i=0; i<obj->pointcnt; i++, ip++) {
        // check if current point is outside of bounding box
        int inside = 1;
		if (ip->x < bbox[0][0] || ip->x > bbox[1][0] ||
			ip->y < bbox[0][1] || ip->y > bbox[1][1])
			inside = 0;
        // if previous point was inside
        if (!out_start) {
            // and current point is inside
            if (inside) {
                // copy current point as reder point
                *op++ = *ip;
                // increase render pointer value
                frames[cwbuf].pnext++;
                if (frames[0].pnext >= frames[0].pmax) { qDebug() << "point overflow"; return; }
            } else {
                // set current point as outside point
                out_start = ip;
                // save current point as last rendered point
                last_render_point = *ip;
			}
        } else if (inside) { // we need to enter from outside point to inside point
			if (!skip_out_start_wait) {
                // render end wait points for last outside point
                for (j=0; j<params.end_wait; j++) {
					*op = *out_start;
					op->color = C_BLACK;
					op++;
					frames[cwbuf].pnext++;
                    if (frames[0].pnext >= frames[0].pmax) { qDebug() << "point overflow"; return; }
                }
			}
			skip_out_start_wait = 0;
            // calculate need tracking distance
            float dx = ip->x - out_start->x;
			float dy = ip->y - out_start->y;
			float distance = fmaxf(fabsf(dx),fabsf(dy));
            // calculate need tracking points number
            int points = ceilf(distance/params.off_speed);
            // if need to track laser
            if (distance > params.snap) {
                // add tracking render points
                for (j=0; j<points; j++) {
					op->x = out_start->x + (dx/(float)points) * j;
					op->y = out_start->y + (dy/(float)points) * j;
					op->color = C_BLACK;
					op++;
					frames[cwbuf].pnext++;
                    if (frames[0].pnext >= frames[0].pmax) { qDebug() << "point overflow"; return; }
                }
                // add start wait render points
                for (j=0; j<params.start_wait; j++) {
					*op = *ip;
					op->color = C_BLACK;
					op++;
					frames[cwbuf].pnext++;
                    if (frames[0].pnext >= frames[0].pmax) { qDebug() << "point overflow"; return; }
                }
            }
            // copy current point as render point
            *op++ = *ip;
            // increase render points pointer value
            frames[cwbuf].pnext++;
            if (frames[0].pnext >= frames[0].pmax) { qDebug() << "point overflow"; return; }
            // reset ouside point
            out_start = NULL;
		}
	}
    // check if last point was inside
    if(!out_start) {
        // duplicate last object point end wait times
        for (i=0; i<params.end_wait; i++) {
			addrndpoint(end->x, end->y, C_BLACK);
		}
        // set last object point as last drawn point
        last_render_point = *end;
	} else {
        // duplicate last point outside end wait times
        for (i=0; i<params.end_wait; i++) {
			addrndpoint(out_start->x, out_start->y, C_BLACK);
		}
        // set last point outside as last drawn point
        last_render_point = *out_start;
	}
}

/// Render frame
float olRenderFrame(int max_fps)
{
	int i;
	int count = 0;

    // calculate minimum render points number for necessary fps
    int min_points = params.rate / max_fps;

    // reset render info
    memset(&last_info, 0, sizeof(last_info));

    // set render points pointer to the first element
    frames[cwbuf].pnext=0;
    // get render objects number
    int cnt = wframe.objcnt;

	float dclosest = 0;
	int clinv = 0;

    // check if object ordering is need
    if (!(params.render_flags & RENDER_NOREORDER)) {
        // render objects according to distance
        ////////////////////////////////////////
        Point closest_to = {-1,-1,0}; // first look for the object nearest the topleft
		while(cnt) {
			Object *closest = NULL;
            // find closest object
            for (i=0; i<wframe.objcnt; i++) {
                // skip "empty" object
                if (!wframe.objects[i].pointcnt)
					continue;
                // skipe "very small" object
                if (wframe.objects[i].pointcnt < params.min_length)
					continue;

                // calculate distance to the object (some strange code, isn't it?)
                float dx = wframe.objects[i].points[0].x - closest_to.x;
				float dy = wframe.objects[i].points[0].y - closest_to.y;
				dx = wframe.objects[i].points[0].x + 1;
				dy = wframe.objects[i].points[0].y + 1;
				float distance = fmaxf(fabsf(dx),fabsf(dy)) + 0.01*(fabsf(dx)+fabsf(dy));
                // check if closest object still was not found or distance is less than previous
                if (!closest || distance < dclosest) {
                    // save closest object
                    closest = &wframe.objects[i];
                    // reset inverse flag
                    clinv = 0;
                    // save distance
                    dclosest = distance;
				}
                // process reverse mode
                if (!(params.render_flags & RENDER_NOREVERSE)) {
                    // calculate distance from the last point
                    dx = wframe.objects[i].points[wframe.objects[i].pointcnt-1].x - closest_to.x;
					dy = wframe.objects[i].points[wframe.objects[i].pointcnt-1].y - closest_to.y;
					distance = fmaxf(fabsf(dx),fabsf(dy)) + 0.01*(fabsf(dx)+fabsf(dy));
                    // check if cloases object stell was not found or distance is less than previous
                    if(!closest || distance < dclosest) {
                        // save closest object
                        closest = &wframe.objects[i];
                        // set inverse flag
                        clinv = 1;
                        // save distance
                        dclosest = distance;
					}
				}
			}
            // check if closest object was found
            if (!closest)
				break;
            // if inverse
            if (clinv) {
                // invert points in the closest object
                Point *pt = closest->points;
				int cnt = closest->pointcnt;
				for (i=0; i<cnt/2; i++) {
					Point tmp = pt[i];
					pt[i] = pt[cnt-i-1];
					pt[cnt-i-1] = tmp;
				}
			}
            // render closest object
            //olLog("%d (%d) (nearest to %f,%f)\n", closest - wframe.objects, closest->pointcnt, closest_to.x, closest_to.y);
			render_object(closest);
			//olLog("[%d] ", frames[cwbuf].pnext);
			//olLog("[LRP:%f %f]\n", last_render_point.x, last_render_point.y);
            // make closest object as "empty"
            closest->pointcnt = 0;
			closest_to = last_render_point;
			cnt--;
            // increase rendered objects count
            last_info.objects++;
		}
		//olLog("\n");
	} else {
        // render all objects
        for (i=0; i<wframe.objcnt; i++) {
            // skip "empty" object
            if (wframe.objects[i].pointcnt < params.min_length)
				continue;
            // render object
            render_object(&wframe.objects[i]);
		}
	}
	wframe.psnext = 0;
	wframe.objcnt = 0;
	count = frames[cwbuf].pnext;
	last_info.points = count;

    // if rendered points number is more than maxinum lenght - need to resample!
    if (params.max_framelen && count > params.max_framelen)
	{
        // source rendered points number
        int in_count = count;
        // target rendered points number
        int out_count = params.max_framelen;
        // check if we have necessary space for resampling
        chkpts(count);

        // pointer to the first point
        Point *pin = frames[cwbuf].points;
        // pointer to the first point for output
        Point *pout = &pin[in_count];

		float pos = 0;
        // scale factor
        float delta = count / (float)out_count;

		count = 0;
		while (pos < (in_count - 1)) {
			int ipos = pos;
			float rest = pos - ipos;

            // interpolate point
            pout->x = pin[ipos].x * (1-rest) + pin[ipos+1].x * rest;
			pout->y = pin[ipos].y * (1-rest) + pin[ipos+1].y * rest;

			if (pin[ipos].color == C_BLACK || pin[ipos+1].color == C_BLACK) {
                // add black point
                pout->color = C_BLACK;
                // move to next point
                pos += 1;
				last_info.resampled_blacks++;
			} else {
                // add color point
                pout->color = pin[ipos].color;
                // increase position by interpolation delta
                pos += delta;
			}

			pout++;
			count++;
		}

		memcpy(pin, &pin[in_count], count * sizeof(*pin));
		frames[cwbuf].pnext = count;
		chkpts(0);
		last_info.resampled_points = count;
	}

    // save last rendered point
    float last_x, last_y;
	if (count) {
		last_x = frames[cwbuf].points[count-1].x;
		last_y = frames[cwbuf].points[count-1].y;
	} else {
		last_x = last_y = 0;
	}
    // if rendered points number is less than minimum
    while(count < min_points) {
        // render last point
        frames[cwbuf].points[count].x = last_x;
		frames[cwbuf].points[count].y = last_y;
		frames[cwbuf].points[count].color = C_BLACK;
		count++;
		last_info.padding_points++;
	}
    // corrent rendered points index
    frames[cwbuf].pnext = count;

    // if audio callback is set
    if (audiocb) {
        // send audio data to the callback
        audiocb(frames[cwbuf].audio_l, frames[cwbuf].audio_r, count);
	} else {
        // reset audio data
        memset(frames[cwbuf].audio_l, 0, sizeof(float)*count);
		memset(frames[cwbuf].audio_r, 0, sizeof(float)*count);
	}

	//olLog("Rendered frame! %d\n", cwbuf);
//	cwbuf = (cwbuf + 1) % fbufs;

    return 0;
}


/////////////////////////////////////////////////////////////////////////
//                         Graphic Operations II
/////////////////////////////////////////////////////////////////////////

/// Set current 2D matrix as identity
void olLoadIdentity(void)
{
    // define static identity matrix
    static const float identity[3][3] = {
		{1,0,0},
		{0,1,0},
		{0,0,1}
	};
    // copy static identity matrix to the current 2D matrix
    memcpy(&mtx2d[0][0], &identity[0][0], sizeof(mtx2d));
}

/// Rotate current 2D matrix at given angle
void olRotate(float theta)
{
    // create rotate matrix
    float rot[9] = {
		cosf(theta),-sinf(theta),0,
		sinf(theta),cosf(theta),0,
		0,0,1,
	};
    // multiply current 2D matrix to rotate matrix
    olMultMatrix(rot);
}

/// Translate current 2D matrix at given offset
void olTranslate(float x, float y)
{
    // create translate matrix
    float trans[9] = {
		1,0,0,
		0,1,0,
		x,y,1,
	};
    // multiply current 2D matrix to translate matrix
    olMultMatrix(trans);
}


/// Scale current 2D matrix at given scale
void olScale(float sx, float sy)
{
    // create scale matrix
    float scale[9] = {
		sx,0,0,
		0,sy,0,
		0,0,1,
	};
    // multiply current 2D matrix to scale matrix
    olMultMatrix(scale);
}

/// Multiply current 2D matrix to the given
void olMultMatrix(float m[9])
{
    // create temporary matrix
    float _new[3][3];

    // multiply current matrix and given
    _new[0][0] = mtx2d[0][0]*m[0] + mtx2d[0][1]*m[1] + mtx2d[0][2]*m[2];
    _new[0][1] = mtx2d[0][0]*m[3] + mtx2d[0][1]*m[4] + mtx2d[0][2]*m[5];
    _new[0][2] = mtx2d[0][0]*m[6] + mtx2d[0][1]*m[7] + mtx2d[0][2]*m[8];
    _new[1][0] = mtx2d[1][0]*m[0] + mtx2d[1][1]*m[1] + mtx2d[1][2]*m[2];
    _new[1][1] = mtx2d[1][0]*m[3] + mtx2d[1][1]*m[4] + mtx2d[1][2]*m[5];
    _new[1][2] = mtx2d[1][0]*m[6] + mtx2d[1][1]*m[7] + mtx2d[1][2]*m[8];
    _new[2][0] = mtx2d[2][0]*m[0] + mtx2d[2][1]*m[1] + mtx2d[2][2]*m[2];
    _new[2][1] = mtx2d[2][0]*m[3] + mtx2d[2][1]*m[4] + mtx2d[2][2]*m[5];
    _new[2][2] = mtx2d[2][0]*m[6] + mtx2d[2][1]*m[7] + mtx2d[2][2]*m[8];

    // copy temporary matrix to the current 2D matrix
    memcpy(&mtx2d[0][0], &_new[0][0], sizeof(mtx2d));
}

/// Push current 2D matrix to the stack
void olPushMatrix(void)
{
    // copy current 2D matrix to the stack active slot
    memcpy(&mtx2ds[mtx2dp][0][0], &mtx2d[0][0], sizeof(mtx2d));
    // increase stack active slot index
    mtx2dp++;
}

/// Pop active 2D matrix from the stack
void olPopMatrix(void)
{
    // decrease stack active slot index
    mtx2dp--;
    // copy current 2D matrix from the stack active slot
    memcpy(&mtx2d[0][0], &mtx2ds[mtx2dp][0][0], sizeof(mtx2d));
}

/// Set current 3D matrix as identity
void olLoadIdentity3(void)
{
    // define static edentity 3D matrix
    static const float identity[4][4] = {
		{1,0,0,0},
		{0,1,0,0},
		{0,0,1,0},
		{0,0,0,1},
	};
    // copy current 3D matrix from the static identity matrix
    memcpy(&mtx3d[0][0], &identity[0][0], sizeof(mtx3d));
}

/// Rotate current 3D matrix over X-axis at the given angle
void olRotate3X(float theta)
{
    // create rotate matrix
    float rot[16] = {
		1,0,0,0,
		0,cosf(theta),sinf(theta),0,
		0,-sinf(theta),cosf(theta),0,
		0,0,0,1
	};
    // multiply current 3D matrix to the rotate matrix
    olMultMatrix3(rot);
}

/// Rotate current 3D matrix over Y-axis at the given angle
void olRotate3Y(float theta)
{
    // create rotate matrix
    float rot[16] = {
		cosf(theta),0,-sinf(theta),0,
		0,1,0,0,
		sinf(theta),0,cosf(theta),0,
		0,0,0,1
	};
    // multiply current 3D matrix to the rotate matrix
    olMultMatrix3(rot);
}

/// Rotate current 3D matrix over Z-axis at the given angle
void olRotate3Z(float theta)
{
    // create rotate matrix
    float rot[16] = {
		cosf(theta), sinf(theta), 0, 0,
		-sinf(theta), cosf(theta), 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1
	};
    // multiply current 3D matrix to the rotate matrix
    olMultMatrix3(rot);
}

/// Translate current 3D matrix at the given offset
void olTranslate3(float x, float y, float z)
{
    // create translate matrix
    float trans[16] = {
		1,0,0,0,
		0,1,0,0,
		0,0,1,0,
		x,y,z,1,
	};
    // multiply current 3D matrix to the translate matrix
    olMultMatrix3(trans);
}


/// Scale current 3D matrix at the given scale
void olScale3(float sx, float sy, float sz)
{
    // create scale matrix
    float trans[16] = {
		sx,0,0,0,
		0,sy,0,0,
		0,0,sz,0,
		0,0,0,1,
	};
    // multiply current 3D matrix to the scale matrix
    olMultMatrix3(trans);
}

/// Multiply current 3D matrix to given
void olMultMatrix3(float m[16])
{
    // create temporary matrix
    float _new[4][4];
    // multiply
    _new[0][0] = mtx3d[0][0]*m[ 0] + mtx3d[0][1]*m[ 1] + mtx3d[0][2]*m[ 2] + mtx3d[0][3]*m[ 3];
    _new[0][1] = mtx3d[0][0]*m[ 4] + mtx3d[0][1]*m[ 5] + mtx3d[0][2]*m[ 6] + mtx3d[0][3]*m[ 7];
    _new[0][2] = mtx3d[0][0]*m[ 8] + mtx3d[0][1]*m[ 9] + mtx3d[0][2]*m[10] + mtx3d[0][3]*m[11];
    _new[0][3] = mtx3d[0][0]*m[12] + mtx3d[0][1]*m[13] + mtx3d[0][2]*m[14] + mtx3d[0][3]*m[15];
    _new[1][0] = mtx3d[1][0]*m[ 0] + mtx3d[1][1]*m[ 1] + mtx3d[1][2]*m[ 2] + mtx3d[1][3]*m[ 3];
    _new[1][1] = mtx3d[1][0]*m[ 4] + mtx3d[1][1]*m[ 5] + mtx3d[1][2]*m[ 6] + mtx3d[1][3]*m[ 7];
    _new[1][2] = mtx3d[1][0]*m[ 8] + mtx3d[1][1]*m[ 9] + mtx3d[1][2]*m[10] + mtx3d[1][3]*m[11];
    _new[1][3] = mtx3d[1][0]*m[12] + mtx3d[1][1]*m[13] + mtx3d[1][2]*m[14] + mtx3d[1][3]*m[15];
    _new[2][0] = mtx3d[2][0]*m[ 0] + mtx3d[2][1]*m[ 1] + mtx3d[2][2]*m[ 2] + mtx3d[2][3]*m[ 3];
    _new[2][1] = mtx3d[2][0]*m[ 4] + mtx3d[2][1]*m[ 5] + mtx3d[2][2]*m[ 6] + mtx3d[2][3]*m[ 7];
    _new[2][2] = mtx3d[2][0]*m[ 8] + mtx3d[2][1]*m[ 9] + mtx3d[2][2]*m[10] + mtx3d[2][3]*m[11];
    _new[2][3] = mtx3d[2][0]*m[12] + mtx3d[2][1]*m[13] + mtx3d[2][2]*m[14] + mtx3d[2][3]*m[15];
    _new[3][0] = mtx3d[3][0]*m[ 0] + mtx3d[3][1]*m[ 1] + mtx3d[3][2]*m[ 2] + mtx3d[3][3]*m[ 3];
    _new[3][1] = mtx3d[3][0]*m[ 4] + mtx3d[3][1]*m[ 5] + mtx3d[3][2]*m[ 6] + mtx3d[3][3]*m[ 7];
    _new[3][2] = mtx3d[3][0]*m[ 8] + mtx3d[3][1]*m[ 9] + mtx3d[3][2]*m[10] + mtx3d[3][3]*m[11];
    _new[3][3] = mtx3d[3][0]*m[12] + mtx3d[3][1]*m[13] + mtx3d[3][2]*m[14] + mtx3d[3][3]*m[15];
    // copy temporary matrix to the current
    memcpy(&mtx3d[0][0], &_new[0][0], sizeof(mtx3d));
}

/// Push current 3D matrix to the stack
void olPushMatrix3(void)
{
    // save current 3D matrix to the stack in active slot
    memcpy(&mtx3ds[mtx3dp][0][0], &mtx3d[0][0], sizeof(mtx3d));
    // increase active slot index
    mtx3dp++;
}

/// Pop current matrix from the stack
void olPopMatrix3(void)
{
    // decrease active slot index
    mtx3dp--;
    // restore current 3D matrix from the active slot
    memcpy(&mtx3d[0][0], &mtx3ds[mtx3dp][0][0], sizeof(mtx3d));
}

/// Transform 3D vector by active 3D matrix
void olTransformVertex3(float *x, float *y, float *z)
{
    // calculate new coodinates
    float px;
	float py;
	float pz;
	float pw;

	px = mtx3d[0][0]**x + mtx3d[0][1]**y + mtx3d[0][2]**z + mtx3d[0][3];
	py = mtx3d[1][0]**x + mtx3d[1][1]**y + mtx3d[1][2]**z + mtx3d[1][3];
	pz = mtx3d[2][0]**x + mtx3d[2][1]**y + mtx3d[2][2]**z + mtx3d[2][3];
	pw = mtx3d[3][0]**x + mtx3d[3][1]**y + mtx3d[3][2]**z + mtx3d[3][3];

    // normalize coordinates
    if (pw == 0) px = py = pz = 0; else {
        px /= pw;
        py /= pw;
        pz /= pw;
    }

    // copy new coordinates
    *x = px;
	*y = py;
	*z = pz;
}

/// Draw 3D vertex
void olVertex3(float x, float y, float z, uint32_t color)
{
    // check args
    if (x != x || y != y || z != z)
    {
        //qDebug() << " x " << x << " y " << y;
        return;
    }

    // check 3D vertex shader
    if(v3shader)
        v3shader(&x, &y, &z, &color); // apply 3D vertex shader

    // transform vector by current 3D matrix
    olTransformVertex3(&x, &y, &z);

    // check again
    if (x != x || y != y)
    {
        //qDebug() << " x " << x << " y " << y;
        return;
    }

    // draw 3D vector
    olVertex(x, y, color);
}

/// Draw rectangle
void olRect(float x1, float y1, float x2, float y2, uint32_t color)
{
    // set drawing mode as LINESTRIP and start drawing
    olBegin(OL_LINESTRIP);
    // draw points
    olVertex(x1,y1,color);
	olVertex(x1,y2,color);
	olVertex(x2,y2,color);
	olVertex(x2,y1,color);
	olVertex(x1,y1,color);
    // stop drawing
    olEnd();
}

/// Draw line
void olLine(float x1, float y1, float x2, float y2, uint32_t color)
{
    // set drawing mode as LINESTRIP and start drawing
    olBegin(OL_LINESTRIP);
    // draw points
    olVertex(x1,y1,color);
	olVertex(x2,y2,color);
    // stop drawing
    olEnd();
}


/// Draw dot
void olDot(float x, float y, int samples, uint32_t color)
{
    // set drawing mode as POINTS and start drawing
    int i;
	olBegin(OL_POINTS);
    // draw point necessary number times
    for (i = 0; i < samples; i++)
		olVertex(x,y,color);
    // stop drawing
    olEnd();
}

/// Reset current color
void olResetColor(void)
{
    // set currnet color as WHITE
    curcol = C_WHITE;
}

/// Multiply current color to the given
void olMultColor(uint32_t color)
{
    // set active color as result of multiplying
    curcol = colmul(curcol, color);
}

/// Push current color to the stack
void olPushColor(void)
{
    // save current color to the stack active slot
    cols[coldp] = curcol;
    // increase stack active slot index
    coldp++;
}

/// Pop current color from the stack
void olPopColor(void)
{
    // decrease stack active slot index
    coldp--;
    // restore current color from the stack active slot
    curcol = cols[coldp];
}


/// Set vertex pre-shader
void olSetVertexPreShader(ShaderFunc f)
{
	vpreshader = f;
}

/// Set vertex shader
void olSetVertexShader(ShaderFunc f)
{
	vshader = f;
}

/// Set vertex 3D shader
void olSetVertex3Shader(Shader3Func f)
{
	v3shader = f;
}

/// Set pixel shader
void olSetPixelShader(ShaderFunc f)
{
	pshader = f;
}

/// Set audio callback
void olSetAudioCallback(AudioCallbackFunc f)
{
	audiocb = f;
}

/// Set frustum for current 3D matrix
void olFrustum (float l, float r, float b, float t, float n, float f)
{
    // create frustum matrix
    float m[16] = {
		(2*n)/(r-l),  0,            0,               0,
		0,            (2*n)/(t-b),  0,               0,
		(r+l)/(r-l),  (t+b)/(t-b),  -(f+n)/(f-n),   -1,
		0,            0,            (-2*f*n)/(f-n),  0,
	};

    // multiply current matrix to the frustum matrix
    olMultMatrix3(m);
}

/// Set perspective parameters for current 3D matrix
void olPerspective (float fovy, float aspect, float zNear, float zFar)
{
    float xmin, xmax, ymin, ymax;

    // calculate maximum value for Y-axis
    ymax = zNear * tanf(fovy * M_PI / 360.0);
    // calculate minimum value for Y-axis
    ymin = -ymax;

    // calculate minimum value for X-axis
    xmin = ymin * aspect;
    // calculate maximum value for X-axis
    xmax = ymax * aspect;

    // set calculated frustum for current 3D matrix
    olFrustum( xmin, xmax, ymin, ymax, zNear, zFar );
}

/// Set drawing bounding box
void olSetScissor (float x0, float y0, float x1, float y1)
{
	bbox[0][0] = x0;
	bbox[0][1] = y0;
	bbox[1][0] = x1;
	bbox[1][1] = y1;
}

/// Get current frame info
void olGetFrameInfo(OLFrameInfo *info)
{
	*info = last_info;
}

/// Log message
void olLog(const char *fmt, ...)
{
    // create buffer
    char buf[1024];
	va_list ap;

    // format message to the buffer
    va_start(ap, fmt);
	vsnprintf(buf, 1024, fmt, ap);
	va_end(ap);
	buf[1023] = 0;

    // check if log callback is set
    if (log_cb)
        // send message to the log callback
        log_cb(buf);
	else
        // send message to the default console
        printf("%s", buf);
}

/// Set log callback
void olSetLogCallback(LogCallbackFunc f)
{
	log_cb = f;
}

RenderedFrame * olGetRenderedFrames(){
    return frames;
}
