#ifndef RPMT_H_
#define RPMT_H_

#include <TCut.h>
const Double_t rpmt_LLD =  200;
const Double_t rpmt_HLD = 7400;
const TCut cut_rpmt_f   = "f==4";
const TCut cut_rpmt_LLD = Form("a>%f && b>%f && c>%f && d>%f",
			       rpmt_LLD, rpmt_LLD, rpmt_LLD, rpmt_LLD);
const TCut cut_rpmt_HLD = Form("a<%f && b<%f && c<%f && d<%f",
			       rpmt_HLD, rpmt_HLD, rpmt_HLD, rpmt_HLD);
const TCut cut_rpmt = cut_rpmt_f && cut_rpmt_LLD && cut_rpmt_HLD;

#endif
