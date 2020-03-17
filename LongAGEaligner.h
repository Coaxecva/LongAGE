/****************************************************************************
 *
 *     LongAGE -- Ultra-long Alignment with Gap Excision
 *     Copyright (C)  Alexej Abyzov
 *                
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Creative Commons license
 * (Attribution-NonCommerical).
 * See terms at http://creativecommons.org/licenses/by-nc/2.5/legalcode
 *
 * Author: Alexej Abyzov
 */

#ifndef __LONG_AGE_ALIGNER_H__
#define __LONG_AGE_ALIGNER_H__

//--- C/C++ includes ---
#include <cstring>
#include <cstdlib>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <limits.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#ifdef LONG_AGE_TIME
#include <sys/time.h>
#endif
using namespace std;

//--- LONG AGE includes ---
#include "Sequence.h"
#include "Scorer.h"
#include "LongAGE.h"

const static unsigned char DIAGONAL   = 1;
const static unsigned char HORIZONTAL = 2;
const static unsigned char VERTICAL   = 3;
const static unsigned char MASK       = 3;

const static unsigned char FS_DIAGONAL   = DIAGONAL;
const static unsigned char FS_HORIZONTAL = HORIZONTAL;
const static unsigned char FS_VERTICAL   = VERTICAL;
const static unsigned char FS_MASK       = MASK;

const static unsigned char RS_DIAGONAL   = DIAGONAL<<2;
const static unsigned char RS_HORIZONTAL = HORIZONTAL<<2;
const static unsigned char RS_VERTICAL   = VERTICAL<<2;
const static unsigned char RS_MASK       = MASK<<2;

const static unsigned char FM_DIAGONAL   = DIAGONAL<<4;
const static unsigned char FM_HORIZONTAL = HORIZONTAL<<4;
const static unsigned char FM_VERTICAL   = VERTICAL<<4;
const static unsigned char FM_MASK       = MASK<<4;

const static unsigned char RM_DIAGONAL   = DIAGONAL<<6;
const static unsigned char RM_HORIZONTAL = HORIZONTAL<<6;
const static unsigned char RM_VERTICAL   = VERTICAL<<6;
const static unsigned char RM_MASK       = MASK<<6;

class LongAGEaligner {
public:
    const static int INDEL_FLAG        = 0x01;
    const static int TDUPLICATION_FLAG = 0x02;
    const static int INVERSION_FLAG    = 0x04;
    const static int INVR_FLAG         = 0x08;
    const static int INVL_FLAG         = 0x10;
    const static int SHOW_ALL_POS      = 0x20;

private:
    const static int ALL_MODES = INDEL_FLAG | TDUPLICATION_FLAG |
                                 INVERSION_FLAG | INVR_FLAG | INVL_FLAG;
    const static int ALL_DISP_OPTS     = SHOW_ALL_POS;

private:
    Sequence &_s1, *_s1_rc, &_s2;
    string &_seq1, &_seq1_rc, &_seq2;
    int _len1, _len2;
    AliFragment *_frags, *_frags_alt;
    int score_n1, score_n2;
    long score_size;
    int _flag;
    int _disp_opts;
    LongAGEaligner *_aux_aligner;
    static const int MAX_BPOINTS = 100;
    int _max;
    short _match, _mismatch, _gap_open, _gap_extend;
public:
    LongAGEaligner(Sequence &s1, Sequence &s2);
    ~LongAGEaligner();
    bool align(Scorer &scr, int flag);
    void printAlignment();
    int score();
private:
    void printMatrix(unsigned short *matr);
    int calcWidth(int v1, int v2);
    int calcOutsideIdentity(string &seq, int bp1, int bp2);
    int calcInsideIdentity(string &seq, int bp1, int bp2);
    int calcIdentity(string &seq, int bp1, int bp2, int &left, int &right);
    void freeFragments(AliFragment *frags);
    bool getExcisedRange(AliFragment *f, int &s, int &e, int add_s = 0, int add_e = 0);
};

#endif