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

#ifndef __LONG_AGE_H__
#define __LONG_AGE_H__

//--- C/C++ includes ---
#include <vector>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#ifdef LONG_AGE_TIME
#include <sys/time.h>
#endif
using namespace std;

#include "Sequence.h"
#include "Scorer.h"

class Location {
public:
    long n;
    long m;
    Location() {n = m = -1;};
    Location(long nn, long mm) {n = nn; m = mm;}
    ~Location() {};
};

class ScorePosition {
public:
    int score;
    long position;
    ScorePosition() {position = score = -1;};
    ~ScorePosition() {};
};

class LongAGEParameters {
public:
    string seq1;
    string seq2;
    int start1, start2;
    int inc1, inc2;
    Scorer* scr;
    bool tdup_inv;
    LongAGEParameters();
    ~LongAGEParameters() {};
};

class AliFragment
{
private:
  AliFragment *_next, *_prev;
  string       _ali1,  _ali2;
  int          _start1,_start2;
  int          _end1,  _end2;
public:
  AliFragment(string &ali1,string &ali2,
          int start1,  int start2,
          int end1,    int end2) : _next(NULL),     _prev(NULL),
                       _ali1(ali1),     _ali2(ali2),
                       _start1(start1), _start2(start2),
                       _end1(end1),     _end2(end2)
  {}

  AliFragment *next(AliFragment *b)
  {
    if (_next) return NULL;
    _next = b;
    b->prev(this);
    return _next;
  }

  AliFragment *prev(AliFragment *b)
  {
    if (_prev) return NULL;
    _prev = b;
    b->next(this);
    return _prev;
  }

  inline AliFragment *next()   { return _next; }
  inline AliFragment *prev()   { return _prev; }
  inline int          start1() { return _start1; }
  inline int          start2() { return _start2; }
  inline int          end1()   { return _end1; }
  inline int          end2()   { return _end2; }

  void countAligned(int &n_ali, int &n_id, int &n_gap)
  {
    n_ali = n_id = n_gap = 0;
    int len = _ali1.length();
    if (_ali2.length() < len) len = _ali2.length();
    for (int i = 0;i < len;i++) {
      n_ali++;
      if (Sequence::sameNuc(_ali1[i],_ali2[i]) > 0)               n_id++;
      if (Sequence::isGap(_ali1[i]) || Sequence::isGap(_ali2[i])) n_gap++;
    }
  }

  void printAlignment()
  {
    static const int WIDTH  = 60;
    static const int MARGIN =  9;

    int inc1 = 1; if (_end1 < _start1) inc1 = -1;
    int inc2 = 1; if (_end2 < _start2) inc2 = -1;

    string margin = "";
    for (int j = 0;j <= MARGIN;j++) margin += ' ';

    int n = _ali1.length(); if (_ali2.length() < n) n = _ali2.length();
    int ind1 = _start1, ind2 = _start2;
    for (int i = 0;i < n;i += WIDTH) {
      int st1 = ind1,st2 = ind2;
      string a1 = _ali1.substr(i,WIDTH);
      string a2 = _ali2.substr(i,WIDTH);
      int nuc1 = 0,nuc2 = 0;
      string match = "";
      for (int j = 0;j < a1.length();j++) {
    if (Sequence::sameNuc(a1[j],a2[j]) > 0) match += '|';
    else if (!Sequence::isGap(a1[j]) && 
         !Sequence::isGap(a2[j]))       match += '.';
    else match += ' ';
    if (!Sequence::isGap(a1[j])) { nuc1++; ind1 += inc1; }
    if (!Sequence::isGap(a2[j])) { nuc2++; ind2 += inc2; }
      }

      cout<<endl;

      if (nuc1 > 0) {
    cout<<setw(MARGIN)<<st1<<' '<<a1;
    cout<<' '<<ind1 - inc1<<endl;
      } else cout<<margin<<a1<<endl;

      cout<<margin<<match<<endl;

      if (nuc2 > 0) {
    cout<<setw(MARGIN)<<st2<<' '<<a2;
    cout<<' '<<ind2 - inc2<<endl;
      } else cout<<margin<<a2<<endl;
    }
  }
};

class LongAGE {
private:
    long score_n2;
    string seq1, seq2;
    bool swap;
    long len1, len2, _score;
    Scorer* scr;
    AliFragment* _frags;
    LongAGEParameters paras;
    Location split(long len_n_m);
    ScorePosition get_fisrt_BP(ScorePosition end);
    pair<string, string> get_aligned_seqs(string s1, string s2);
    pair<ScorePosition, pair<ScorePosition, long> > local_run(long n_stop, long m_start, long m_stop, int* M);

public:
    long Score();
    AliFragment* Fragments();
    LongAGE(LongAGEParameters paras);
    void findAlignment(LongAGEParameters paras);
    pair<pair<long, long>,pair<long,long> > find_BPs(Scorer* scr, bool split_matrix);
    pair<pair<long, long>,pair<long,long>> find_inv_tdup_BPs(Scorer* scr);
    pair<pair<long, long>,pair<long,long>> find_indel_BPs();
    pair<pair<ScorePosition, ScorePosition>,pair<ScorePosition,ScorePosition> > get_candidatate_BPs(Scorer* scr, int* M);
    ~LongAGE();
};

#endif