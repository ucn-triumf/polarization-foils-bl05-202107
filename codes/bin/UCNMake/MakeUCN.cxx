#include "MakeUCN.h"
#include "DrawNikiRPMT.cxx"

extern void InitGui();

VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
TROOT api("nikiglass","nikiglass display",initfuncs);

int main(int argc, char *argv[])
{

	gStyle->SetPalette(1);
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(1001111);//to delete number pallete
	gStyle->SetOptFit(1110);
	TGaxis::SetMaxDigits(3);

	// ******************************************************************************
	// timer for information about running time of code
	TStopwatch timer;
	timer.Start();

	std::ostringstream os, osopt;
	char DataFileName[128];

	if(argc==1 || argc>3){
		std::cerr<<"Usage: ./MakeNiki [DataFileName] [Option Make Both:0(default), Only RPMT:1, Only single:2]" <<std::endl;
		return -1;
	} else {
		os <<argv[1];
		sprintf(DataFileName,os.str().c_str());
		std::cout <<"Making ROOT file from "<<DataFileName<<std::endl;
		if(argc==2){
			nOption = 0;
		}else if(argc==3){
			osopt <<argv[2];
			nOption = atoi(osopt.str().c_str());
			std::cout <<"Option " <<nOption <<" was chosen."<<std::endl;
		} 
	}

	if(nOption==0)   std::cout <<"Create RPMT and Single channel ntuples."<<std::endl;
	if(nOption==1)   std::cout <<"Create Only a RPMT ntuples."<<std::endl;
	if(nOption==2)   std::cout <<"Create Only a Single channel ntuple."<<std::endl;

	MakeNikiRPMT(DataFileName); //for RPMT 
	if(nOption!=2)  DrawNikiRPMT(DataFileName); //for RPMT 

	//  MakeIwatsuCh(DataFileName); //for beam monitor

	// ******************************************************************************
	// timer for information about running time of code
	cout<<"*********************************"<<endl;
	cout<<"ROOT - Time at the end of job:"<<endl;
	cout<<"CpuTime = "<<timer.CpuTime()<<" seconds"<<endl;
	cout<<"RealTime = "<<timer.RealTime()<<" seconds"<<endl;
	cout<<"*********************************"<<endl<<flush;
	timer.Stop();

	return 0;

}

//for single channel data
Int_t *MakeNikiCh(char* DataFileName)
{
	Int_t nlines=0;
	Int_t kp=0,tp=0,timing=0;
	Int_t board,ch,e,dummy; 
	long long int t=0,tzero=0;    // use long long t t!!
	long long int Tzero=0;    // use long long t t!!
	TNtuple *seven = new  TNtuple("Seven","Seven","e:t:tof:timing:kp::tp:ch:board"); //"f" is coincidence number of event: use f==4

	FILE *fp;
	fp = fopen(DataFileName, "r");
	if (!fp) {
		cerr << "File could not be opened \n"<<DataFileName<<endl;
		return 0;
	}

	// data read start
	while(fscanf(fp, "%d,%d,%d,%lld,%d",&board,&ch,&e,&t,&dummy)!=EOF){

		if (nlines%1000000==0) cout <<" nlines =  "<<nlines<<endl;
		//for noise reduction
		if(1){
			//for kicker pulse (beam timing)
			//if data is kp trigger, then reset tzero to 0.
			if(ch==timing_ch && e > kp_LLD && e<kp_HLD){
				tzero=t;
				timing++;
#if KPROOT
				seven->Fill(e,t,t-tzero,timing,kp,tp,ch,board);
#endif
			}
			else if(ch==kp_ch && e > kp_LLD && e<kp_HLD){
				Tzero=t;
				kp++;
#if KPROOT
				//add kp in root
				seven->Fill(e,t,t-tzero,timing,kp,tp,ch,board);
#endif
			}
			//for timing pulse(25Hz)
			else if(ch==tp_ch && e > kp_LLD && e<kp_HLD){
				//	tzero=t;
				tp++;
#if KPROOT
				//add tp in root
				seven->Fill(e,t,t-tzero,timing,kp,tp,ch,board);
#endif
			}
			else{
				seven->Fill(e,t,t-tzero,timing,kp,tp,ch,board);
			}
		}
		nlines++;
	}

	//for last event
	if(ch==timing_ch && e > kp_LLD && e<kp_HLD){
		timing++;
#if KPROOT
		seven->Fill(e,t,t-tzero,timing,kp,tp,ch,board);
#endif
	}else if(ch==kp_ch){
		kp++;
#if KPROOT
		seven->Fill(e,t,t-tzero,timing,kp,tp,ch,board);
#endif
	}
	else if(ch==tp_ch){
		tp++;
#if KPROOT
		seven->Fill(e,t,t-tzero,timing,kp,tp,ch,board);
#endif
	}
	else seven->Fill(e,t,t-tzero,timing,kp,tp,ch,board);
	fclose(fp);
	cout <<"KP= "<<kp<<endl;
	cout <<"Timing Pilse From Doppler Shifter = "<<timing<<endl;

	//create tof.txt file for Taketani automatic scan system.
	FILE *tofp;
	tofp = fopen("/home/iwatsu/bin/tof.txt", "w");
	if (!tofp) {
		cerr << "/home/iwatsu/bin/tof.txt could not be opened"<<endl;
		return 0;
	}

	fprintf(tofp,"%d\n",kp); // number of T0 counts
	fprintf(tofp,"100\n"); // bin width is 100 usec

	TH1F *htof = new TH1F("htof","TOF",400,0,40);
	seven->Draw("tof*1e-3>>htof",Form("e > %d && ch == %d", monitor_thr, monitor_ch),"off");
	int n = htof->GetNbinsX();
	float *tmp = htof->GetArray();
	for(int i=0; i<n; i++){
		fprintf(tofp,"%lf\n",tmp[i]); 
	}
	fclose(tofp);

	char RootFileName[64];
	sprintf(RootFileName,"%s.kp.root",DataFileName);
	TFile file(RootFileName,"RECREATE");
	seven->Write();
	cout << RootFileName << " has been created" << endl;
	file.Close();

cout << "aaaaaaaaaaaaaaaaaaaaaaaa" << endl;

	return 0;
}


//for RPMT and other data
Int_t *MakeNikiRPMT(char* DataFileName)
{
	Int_t i=0,nlines=0;
	Int_t board,ch,e,dummy;  
	Int_t kp=0,tp=0,timing=0;
	Int_t t0flag=0,t1flag=0;
	long long int t=0;            // use long long t t!!
	long long int tzero=0;        // use long long t t!!
	long long int Tzero=0;        // use long long t t!!
	Double_t dt=1.1; //coincidence time for RPMT channels

	struct sCntl sRpmt;
	sRpmt=ClearCntl(sRpmt);

	//RPMT root file
	TNtuple *tup = new  TNtuple("T","T","a:b:c:d:t:tof:f:timing:kp:tp:x:y"); //"f" is coincidence number of event: use f==4
	//other root file
	TNtuple *seven = new  TNtuple("Seven","Seven","e:t:tof:timing:kp:tp:ch:board"); //"f" is coincidence number of event: use f==4

	FILE *fp;
	fp = fopen(DataFileName, "r");
	if (!fp) {
		cerr << "File could not be opened \n"<<DataFileName<<endl;
		return 0;
	}

	// data read start
	while(fscanf(fp, "%d,%d,%d,%lld,%d",&board,&ch,&e,&t,&dummy)!=EOF){

		if (nlines%1000000==0){
			cout <<" nlines =  "<<nlines<<", timing pulse from DS = "<<timing<<endl;
		}

		//This counts the number of DS Timing pulse.
		if(ch==timing_ch && e > kp_LLD && e<kp_HLD){
			timing++;
			tzero=t; 
			t0flag=1;
			//This counts the number of timing pulse.
		}else if(ch==tp_ch && e > kp_LLD && e<kp_HLD){
			tp++;
			t1flag=1;
			//This counts the number of timing pulse.
		}else if(ch==kp_ch && e > kp_LLD && e<kp_HLD){
			Tzero=t; 
			kp++;
		}

		//for other ch events

		//if(ch==7)
		//	printf("t=%lld,tzero=%lld\n",t,tzero);

		if(ch==6){
			seven->Fill(e,t,t-tzero,timing,kp,tp,ch,board);
			seven->Fill(e,t,t-Tzero,-1,kp,tp,ch,board);
		}else{
			seven->Fill(e,t,t-tzero,timing,kp,tp,ch,board);
		}

		//ignore events before first kp.
		if(timing==0 && t0flag==0){
			;
		}
		//if time is later than sett+dt, then Fill data in tuple
		else if(t-sRpmt.sett>dt || t0flag==1 || t1flag==1){

			Double_t RPMT_X, RPMT_Y;
			for(i=0;i<4;i++) sRpmt.iSumFlag+=sRpmt.iFlag[i];
			if((sRpmt.iPulse[0]+sRpmt.iPulse[1])!=0)  RPMT_X = Double_t(sRpmt.iPulse[0])/Double_t(sRpmt.iPulse[0]+sRpmt.iPulse[1]);
			else  RPMT_X = -1.;
			if((sRpmt.iPulse[2]+sRpmt.iPulse[3])!=0)  RPMT_Y = Double_t(sRpmt.iPulse[2])/Double_t(sRpmt.iPulse[2]+sRpmt.iPulse[3]);
			else  RPMT_Y = -1.;
			tup->Fill(sRpmt.iPulse[0],sRpmt.iPulse[1],sRpmt.iPulse[2],sRpmt.iPulse[3],
					sRpmt.sett,sRpmt.sett-tzero,sRpmt.iSumFlag,timing,kp,tp,RPMT_X,RPMT_Y);

			t0flag=0;
			t1flag=0;

			// initialize structure
			sRpmt = ClearCntl(sRpmt);

			//reseting sett
			sRpmt.sett=t;
		}

		//fill new data
		Int_t rpmtch = -1;
		if(ch==rpmt_ch1) rpmtch=0;
		else if(ch==rpmt_ch2) rpmtch=1;
		else if(ch==rpmt_ch3) rpmtch=2;
		else if(ch==rpmt_ch4) rpmtch=3;
		if(rpmtch >= 0){
			sRpmt.iFlag[rpmtch]++;
			sRpmt.iPulse[rpmtch]=e;
		}

		//next line
		nlines++;   
	}

	fclose(fp);

	//for last is not taken into account

#define CREATE_ROOT 1
	char RootFileName[64];
#if CREATE_ROOT
	sprintf(RootFileName,"%s.rpmt.root",DataFileName);
	TFile file(RootFileName,"RECREATE");
	if(nOption!=2)tup->Write();
	if(nOption!=1)seven->Write();
	cout << RootFileName << " has been created" << endl;
	file.Close();
#endif

	std::cout << "timing pulse from DS= " << timing << std::endl;


	return 0;
}

