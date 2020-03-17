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

#include <vector>
#include <algorithm>
#include "LongAGEaligner.h"
#include "LongAGE.h"

LongAGEaligner::LongAGEaligner(Sequence &s1,Sequence &s2) :
        _s1(s1), _s2(s2),
        _seq1(_s1.sequence()), _len1(_seq1.length()),
        _seq2(_s2.sequence()), _len2(_seq2.length()),
        _s1_rc(_s1.clone()),_seq1_rc(_s1_rc->sequence()),
        _frags(NULL), _frags_alt(NULL),_max(0), _flag(0), _aux_aligner(NULL),
        _match(0), _mismatch(0), _gap_open(0), _gap_extend(0)
{
    _s1_rc->revcom();
    score_n1 = _len1 + 2;
    score_n2 = _len2 + 2;
}

LongAGEaligner::~LongAGEaligner()
{
    delete _s1_rc;
    freeFragments(_frags);
    freeFragments(_frags_alt);
    delete _aux_aligner;
}

void LongAGEaligner::freeFragments(AliFragment *frags)
{
    AliFragment *f = frags;
    while (f) {
        AliFragment *tmp = f;
        f = f->next();
        delete tmp;
    }
}

bool LongAGEaligner::align(Scorer &scr,int flag)
{
    if (_len1 <= 0 || _len2 <= 0) return false;
    if (!(flag & ALL_MODES)) return false;

    // Deciding on whether need to use auxiliary aligner
    int aux_flag = 0;
    _flag = flag & ALL_MODES;
    _disp_opts = flag & ALL_DISP_OPTS;
    _flag = flag;
    if (flag & INVERSION_FLAG) {
        _flag    = INVL_FLAG;
        aux_flag = INVR_FLAG  | _disp_opts;
    }

    _match      = scr.getMatch();
    _mismatch   = scr.getMismatch();
    _gap_open   = scr.getGapOpen();
    _gap_extend = scr.getGapExtend(); 


    LongAGEParameters paras;
    paras.scr = &scr;
    paras.seq1 = flag & INVERSION_FLAG ? _seq1_rc : _seq1;
    paras.seq2 = _seq2;
    paras.start1 = _s1.start();
    paras.start2 = _s2.start();    
    paras.inc1 = 1; if (_s1.reverse()) paras.inc1 = -1;
    paras.inc2 = 1; if (_s2.reverse()) paras.inc2 = -1;
    paras.tdup_inv = ((bool)(flag & TDUPLICATION_FLAG)) || ((bool)(flag & INVERSION_FLAG));

    LongAGE LongAGE(paras);
    LongAGE.findAlignment(paras);
    _max = LongAGE.Score();
    _frags = LongAGE.Fragments();

    delete _aux_aligner;
    _aux_aligner = NULL;
    if (aux_flag) {
        _aux_aligner = new LongAGEaligner(_s1,_s2);
        if (!_aux_aligner->align(scr,aux_flag)) {
            delete _aux_aligner;
            _aux_aligner = NULL;
        }
    }

    return true;
}

int LongAGEaligner::score()
{
    if (_aux_aligner && _aux_aligner->score() > _max)
        return _aux_aligner->score();
    return _max;
}

bool LongAGEaligner::getExcisedRange(AliFragment *f,int &s,int &e,
                                 int add_s,int add_e)
{
    if (!f->next()) return false;
    int inc1 = 1; if (_s1.reverse()) inc1 = -1;
    if (_flag & INDEL_FLAG) {
        s = f->end1() + inc1;
        e = f->next()->start1() - inc1;
    } else if (_flag & TDUPLICATION_FLAG) {
        s = f->next()->start1();
        e = f->end1();
    } else if (_flag & INVL_FLAG) {
        s = f->end1() + inc1;
        e = f->next()->start1();
    } else if (_flag & INVR_FLAG) {
        s = f->end1();
        e = f->next()->start1() - inc1;
    } else return false;

    s += inc1*add_s;
    e += inc1*add_e;

    return true;
}

void LongAGEaligner::printAlignment()
{
    if (_aux_aligner && _aux_aligner->score() > _max) {
        _aux_aligner->printAlignment();
        return;
    }

    const static string EXCISED_MESSAGE  = "EXCISED REGION";
    const static string EXCISED_MESSAGEs = "EXCISED REGION(S)";

    cout<<endl;
    cout<<"MATCH = "<<_match<<", ";
    cout<<"MISMATCH = "<<_mismatch<<", ";
    cout<<"GAP OPEN = "<<_gap_open<<", ";
    cout<<"GAP EXTEND = "<<_gap_extend;
    if      (_flag & INDEL_FLAG)        cout<<", INDEL";
    else if (_flag & INVL_FLAG)         cout<<", INVERSION";
    else if (_flag & INVR_FLAG)         cout<<", INVERSION";
    else if (_flag & TDUPLICATION_FLAG) cout<<", TDUPLICATION";
    cout<<endl<<endl;

    int inc1 = 1; if (_s1.reverse()) inc1 = -1;
    int inc2 = 1; if (_s2.reverse()) inc2 = -1;
    int s1 = _s1.start(),s2 =_s2.start();

    int e1 = s1 + inc1*(_seq1.length() - 1);
    int e2 = s2 + inc2*(_seq2.length() - 1);
    cout<<"First  seq ["
        <<setw(calcWidth(s1,s2))<<s1<<","
        <<setw(calcWidth(e1,e2))<<e1<<"] => "
        <<setw(9)<<_seq1.length()<<" nucs '"<<_s1.name()<<"'"<<endl;
    cout<<"Second seq ["
        <<setw(calcWidth(s1,s2))<<s2<<","
        <<setw(calcWidth(e1,e2))<<e2<<"] => "
        <<setw(9)<<_seq2.length()<<" nucs '"<<_s2.name()<<"'"<<endl<<endl;

    int n_frg = 0;
    for (AliFragment *f = _frags;f;f = f->next()) n_frg++;
    int *n_ali = new int[n_frg + 1];
    int *n_id  = new int[n_frg + 1];
    int *n_gap = new int[n_frg + 1];
    n_ali[0] = n_id[0] = n_gap[0] = 0;
    int index = 1;
    for (AliFragment *f = _frags;f;f = f->next()) {
        f->countAligned(n_ali[index],n_id[index],n_gap[index]);
        n_ali[0] += n_ali[index];
        n_id[0]  += n_id[index];
        n_gap[0] += n_gap[index];
        index++;
    }

    int identic = 0,gap = 0;
    if (n_ali[0] > 0) {
        identic = (int)(100.*n_id[0]/n_ali[0] + 0.5);
        gap     = (int)(100.*n_gap[0]/n_ali[0] + 0.5);
    }

    cout<<"Score:   "<<setw(9)<<_max<<endl;
    cout<<"Aligned: "<<setw(9)<<n_ali[0]<<"        nucs"<<endl;
    cout<<"Identic: "<<setw(9)<<n_id[0]<<" ("<<setw(3)<<identic<<"%) nucs";
    if (n_frg > 1) {
        cout<<" =>";
        for (int fr = 1;fr < n_frg + 1;fr++)
            cout<<" "<<setw(9)<<n_id[fr]<<" ("<<setw(3)
            <<(int)(100.*n_id[fr]/n_ali[fr] + 0.5)<<"%)";
    }
    cout<<endl;
    cout<<"Gaps:    "<<setw(9)<<n_gap[0]<<" ("<<setw(3)
    <<gap<<"%) nucs"<<endl<<endl;

    // Printing aligned region coordinates
    if (n_ali[0] > 0) {
        cout<<"Alignment:"<<endl;
        cout<<" first  seq =>  ";
        for (AliFragment *f = _frags;f;f = f->next()) {
            if (f->prev()) cout<<" "<<EXCISED_MESSAGE<<" ";
            int ws = calcWidth(f->start1(),f->start2());
            int we = calcWidth(f->end1(),f->end2());
            cout<<"["<<setw(ws)<<f->start1()<<","<<setw(we)<<f->end1()<<"]";
        }
        cout<<endl;
        cout<<" second seq =>  ";
        for (AliFragment *f = _frags;f;f = f->next()) {
            if (f->prev()) cout<<" "<<EXCISED_MESSAGE<<" ";
            int ws = calcWidth(f->start1(),f->start2());
            int we = calcWidth(f->end1(),  f->end2());
            cout<<"["<<setw(ws)<<f->start2()<<","<<setw(we)<<f->end2()<<"]";

        }
        cout<<endl;
    }

    int s,e,len;
    if (n_frg > 1) {
        cout<<endl<<EXCISED_MESSAGEs<<":"<<endl;
        for (AliFragment *f = _frags;f && f->next();f = f->next()) {
            if (!getExcisedRange(f,s,e)) continue;
            cout<<" first  seq => ";
            len = abs(s - e - inc1);
            cout<<setw(9)<<len<<" nucs";
            if (len > 0) cout<<" ["<<s<<","<<e<<"]";
            cout<<endl;
            cout<<" second seq => ";
            s = f->end2() + inc2;
            e = f->next()->start2() - inc2;
            len = abs(s - e - inc2);
            cout<<setw(9)<<len<<" nucs";
            if (len > 0) cout<<" ["<<s<<","<<e<<"]";
            cout<<endl;
        }

        // Printing sequence identity around breakpoints
        cout<<endl<<"Identity at breakpoints: "<<endl;
        int left,right,s,e;

        for (AliFragment *f = _frags;f && f->next();f = f->next()) {
            if (!getExcisedRange(f,s,e,0,1)) continue;
            int bp1 = inc1*(s - _s1.start()), bp2 = inc1*(e - _s1.start());
            int n_hom = calcIdentity(_seq1,bp1,bp2,left,right);
            cout<<" first  seq => "<<setw(9)<<n_hom<<" nucs";
            if (n_hom > 0)
                cout<<" ["<<s + inc1*left<<","<<s + inc1*(right - 1)<<"] to"
                <<" ["<<e + inc1*left<<","<<e + inc1*(right - 1)<<"]";
            cout<<endl;
            s = f->end2() + inc2, e = f->next()->start2();
            bp1 = inc2*(s - _s2.start()), bp2 = inc2*(e - _s2.start());
            n_hom = calcIdentity(_seq2,bp1,bp2,left,right);
            cout<<" second seq => "<<setw(9)<<n_hom<<" nucs";
            if (n_hom > 0)
                cout<<" ["<<s + inc2*left<<","<<s + inc2*(right - 1)<<"] to"
                <<" ["<<e + inc2*left<<","<<e + inc2*(right - 1)<<"]";
            cout<<endl;
        }

        cout<<"Identity outside breakpoints: "<<endl;
        for (AliFragment *f = _frags;f && f->next();f = f->next()) {
            if (!getExcisedRange(f,s,e,-1,1)) continue;
            int bp1 = inc1*(s - _s1.start()), bp2 = inc1*(e - _s1.start());
            int n_hom = calcOutsideIdentity(_seq1,bp1,bp2);
            cout<<" first  seq => "<<setw(9)<<n_hom<<" nucs";
            if (n_hom > 0)
                cout<<" ["<<s - inc1*(n_hom - 1)<<","<<s<<"] to"
                <<" ["<<e<<","<<e + inc1*(n_hom - 1)<<"]";
            cout<<endl;
            s = f->end2(), e = f->next()->start2();
            bp1 = inc2*(s - _s2.start()), bp2 = inc2*(e - _s2.start());
            n_hom = calcOutsideIdentity(_seq2,bp1,bp2);
            cout<<" second seq => "<<setw(9)<<n_hom<<" nucs";
            if (n_hom > 0)
                cout<<" ["<<s - inc2*(n_hom - 1)<<","<<s<<"] to"
                <<" ["<<e<<","<<e + inc2*(n_hom - 1)<<"]";
            cout<<endl;
        }

        cout<<"Identity inside breakpoints: "<<endl;
        for (AliFragment *f = _frags;f && f->next();f = f->next()) {
            if (!getExcisedRange(f,s,e)) continue;
            int bp1 = inc1*(s - _s1.start()), bp2 = inc1*(e - _s1.start());
            int n_hom = calcInsideIdentity(_seq1,bp1,bp2);
            cout<<" first  seq => "<<setw(9)<<n_hom<<" nucs";
            if (n_hom > 0)
                cout<<" ["<<s<<","<<s + inc1*(n_hom - 1)<<"] to"
                <<" ["<<e - inc1*(n_hom - 1)<<","<<e<<"]";
            cout<<endl;
            s = f->end2() + inc2, e = f->next()->start2() - inc2;
            bp1 = inc2*(s - _s2.start()), bp2 = inc2*(e - _s2.start());
            n_hom = calcInsideIdentity(_seq2,bp1,bp2);
            cout<<" second seq => "<<setw(9)<<n_hom<<" nucs";
            if (n_hom > 0)
                cout<<" ["<<s<<","<<s + inc2*(n_hom - 1)<<"] to"
                <<" ["<<e - inc2*(n_hom - 1)<<","<<e<<"]";
            cout<<endl;
        }
    }

    // Printing actual alignment
    for (AliFragment *f = _frags;f;f = f->next()) {
        if (f != _frags) cout<<endl<<EXCISED_MESSAGE<<endl;
        f->printAlignment();
    }

    delete[] n_ali;
    delete[] n_id;
    delete[] n_gap;
}

int LongAGEaligner::calcWidth(int v1,int v2)
{
    int width1 = 0; if (v1 < 0) v1++;
    while (v1 != 0) { width1++; v1 /= 10; }
    int width2 = 0; if (v2 < 0) v2++;
    while (v2 != 0) { width2++; v2 /= 10; }
    if (width1 > width2) return width1;
    return width2;
}

int LongAGEaligner::calcOutsideIdentity(string &seq,int bp1,int bp2)
{
    // bp1 and bp2 are zero based
    if (bp1 >= bp2) return 0;
    int n = seq.length(),d1 = bp1 + 1, d2 = n - bp2;
    int n_check = d1; if (d2 < d1) n_check = d2;
    int start1 = bp1 - n_check + 1;
    for (int i = 0;i < n_check;i++) {
        int n_same = 0,start2 = bp2 - i;
        for (int j = i;j < n_check;j++)
            if (Sequence::sameNuc(seq[start1 + j],seq[start2 + j])) n_same++;
            else break;
        if (n_same == n_check - i) return n_same;
    }

    return 0;
}

int LongAGEaligner::calcInsideIdentity(string &seq,int bp1,int bp2)
{
    // bp1 and bp2 are zero based
    if (bp1 >= bp2) return 0;
    int n_check = (bp2 - bp1 + 1)/2, start2 = bp2 - n_check + 1;
    for (int i = 0;i < n_check;i++) {
        int n_same = 0,start1 = bp1 - i;
        for (int j = i;j < n_check;j++) {
            if (Sequence::sameNuc(seq[start1 + j],seq[start2 + j])) n_same++;
            else break;
        }
        if (n_same == n_check - i) return n_same;
    }

    return 0;
}

int LongAGEaligner::calcIdentity(string &seq,int bp1,int bp2,int &left,int &right)
{
    // bp1 and bp2 are zero based
    if (bp1 >= bp2) return 0;
    int n = seq.length(),delta = bp2 - bp1;
    if (delta == 0) return 0;
    int start = bp1 - 1;
    left = right = 0;
    while (start >= 0 && Sequence::sameNuc(seq[start],seq[start + delta])) {
        start--;
        left--;
    }

    int end = bp1;
    while (end < n && Sequence::sameNuc(seq[end],seq[end + delta])) {
        end++;
        right++;
    }

    return right - left;
}

void LongAGEaligner::printMatrix(unsigned short *matr)
{
    long pos = 0;
    for (int i2 = 0;i2 < score_n2;i2++) {
        for (int i1 = 0;i1 < score_n1;i1++)
            cout<<setw(4)<<matr[pos++]<<" ";
        cout<<endl;
    }
    cout<<endl;
}
