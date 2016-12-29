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

#ifndef LIBOL_H
#define LIBOL_H

#include <stdint.h>

#include <QtCore/qglobal.h>

#ifdef Q_OS_WIN32
#define OPENLASE_EXPORT __declspec(dllexport)
#else
#define OPENLASE_EXPORT
#endif


/////////////////////////////////////////////////////////////////////////
//                         Data Structures
/////////////////////////////////////////////////////////////////////////





/*
    struct for points
    used in struct Object, struct Frame, strcut RenderedFrame
    float x, y, range in [-1, 1]
    uint32_t color, 8bit R,G,B each, top byte is unused
*/
struct Point {
    float x,y;
    uint32_t color;
};

/*
    struct RenderedFrame, finished frame
    int pmax, maximum points allowed for this frame
    int pnext, offset to last+1 point of points array
    Point *points, pointer to array of points
    float *audio_l, left channel audio
    float *audio_r, right chanell audio
*/
struct RenderedFrame {
    int pmax;
    int pnext;
    Point *points;
    float *audio_l;
    float *audio_r;
};

enum {
	OL_LINESTRIP,
	OL_BEZIERSTRIP,
	OL_POINTS
};

#define C_RED   0xff0000
#define C_GREEN 0x00ff00
#define C_BLUE  0x0000ff
#define C_WHITE 0xffffff
#define C_BLACK 0x000000

#define C_GREY(x)   (0x010101 * ((int)(x)))
#define C_RED_I(x)   (0x010000 * ((int)(x)))
#define C_GREEN_I(x)   (0x000100 * ((int)(x)))
#define C_BLUE_I(x)   (0x000001 * ((int)(x)))


enum {
	RENDER_GRAYSCALE = 1,
	RENDER_NOREORDER = 2,
	RENDER_NOREVERSE = 4,
};

struct OLRenderParams{
	int rate;
	float on_speed;
	float off_speed;
	int start_wait;
	int start_dwell;
	int curve_dwell;
	int corner_dwell;
	int end_dwell;
	int end_wait;
	float curve_angle;
	float flatness;
	float snap;
	int render_flags;
	int min_length;
	int max_framelen;
};


struct OLFrameInfo {
	int objects;
	int points;
	int resampled_points;
	int resampled_blacks;
	int padding_points;
};

/////////////////////////////////////////////////////////////////////////
//                         Libol Module Functions
/////////////////////////////////////////////////////////////////////////
OPENLASE_EXPORT int olInit(int max_points);

OPENLASE_EXPORT void olSetRenderParams(OLRenderParams *params);
OPENLASE_EXPORT void olGetRenderParams(OLRenderParams *params);

typedef void (*AudioCallbackFunc)(float *leftbuf, float *rightbuf, int samples);

void olSetAudioCallback(AudioCallbackFunc f);


OPENLASE_EXPORT float olRenderFrame(int max_fps);

OPENLASE_EXPORT void olGetFrameInfo(OLFrameInfo *info);

OPENLASE_EXPORT RenderedFrame * olGetRenderedFrames();

OPENLASE_EXPORT void olSetScissor (float x0, float y0, float x1, float y1);

OPENLASE_EXPORT void olLog(const char *fmt, ...);

typedef void (*LogCallbackFunc)(const char *msg);

OPENLASE_EXPORT void olSetLogCallback(LogCallbackFunc f);
/////////////////////////////////////////////////////////////////////////
//                         Libol Shader Functions
/////////////////////////////////////////////////////////////////////////

typedef void (*ShaderFunc)(float *x, float *y, uint32_t *color);
typedef void (*Shader3Func)(float *x, float *y, float *z, uint32_t *color);

OPENLASE_EXPORT void olSetVertexPreShader(ShaderFunc f);
OPENLASE_EXPORT void olSetVertexShader(ShaderFunc f);
OPENLASE_EXPORT void olSetVertex3Shader(Shader3Func f);

OPENLASE_EXPORT void olSetPixelShader(ShaderFunc f);


/////////////////////////////////////////////////////////////////////////
//                         Libol Graphic Operations
/////////////////////////////////////////////////////////////////////////

OPENLASE_EXPORT void olLoadIdentity(void);
OPENLASE_EXPORT void olPushMatrix(void);
OPENLASE_EXPORT void olPopMatrix(void);

OPENLASE_EXPORT void olMultMatrix(float m[9]);
OPENLASE_EXPORT void olRotate(float theta);
OPENLASE_EXPORT void olTranslate(float x, float y);
OPENLASE_EXPORT void olScale(float sx, float sy);


OPENLASE_EXPORT void olLoadIdentity3(void);
OPENLASE_EXPORT void olPushMatrix3(void);
OPENLASE_EXPORT void olPopMatrix3(void);

OPENLASE_EXPORT void olMultMatrix3(float m[16]);
OPENLASE_EXPORT void olRotate3X(float theta);
OPENLASE_EXPORT void olRotate3Y(float theta);
OPENLASE_EXPORT void olRotate3Z(float theta);
OPENLASE_EXPORT void olTranslate3(float x, float y, float z);
OPENLASE_EXPORT void olScale3(float sx, float sy, float sz);

OPENLASE_EXPORT void olFrustum (float left, float right, float bot, float ttop, float near, float far);
OPENLASE_EXPORT void olPerspective(float fovy, float aspect, float zNear, float zFar);

OPENLASE_EXPORT void olResetColor(void);
OPENLASE_EXPORT void olMultColor(uint32_t color);
OPENLASE_EXPORT void olPushColor(void);
OPENLASE_EXPORT void olPopColor(void);

OPENLASE_EXPORT void olBegin(int prim);
OPENLASE_EXPORT void olVertex(float x, float y, uint32_t color);
OPENLASE_EXPORT void olVertex3(float x, float y, float z, uint32_t color);
OPENLASE_EXPORT void olEnd(void);

OPENLASE_EXPORT void olTransformVertex3(float *x, float *y, float *z);


OPENLASE_EXPORT void olRect(float x1, float y1, float x2, float y2, uint32_t color);
OPENLASE_EXPORT void olLine(float x1, float y1, float x2, float y2, uint32_t color);
OPENLASE_EXPORT void olDot(float x, float y, int points, uint32_t color);



#endif
