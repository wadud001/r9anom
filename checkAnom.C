#include "extra_tools.cc"
#include "TChain.h"
#include <iostream>
#include <fstream>

std::string NTUPLE_LIST = "UL_baseline.txt";
std::string OUT_FILE = "r9_anomalies.root";

std::vector<Short_t> INDICES = {0,1};

void getAlphaStuff(std::string _ntuple_list = NTUPLE_LIST, std::string _out_file = OUT_FILE){
	std::cout<<"Using ntuple list "<<_ntuple_list<<" and will write to "<<_out_file<<std::endl;
	std::vector<std::string> _ntuples = getNonemptyLines(_ntuple_list);
	TChain *_ntuple_chain = new TChain("selected");
	for(auto & ntuple : _ntuples){
		_ntuple_chain->Add(ntuple.c_str());
		std::cout << "\t" << ntuple << std::endl;
	}
	ULong64_t nEntries = _ntuple_chain->GetEntries();
	_ntuple_chain->LoadTree(-1);
	ULong64_t reportEvery = std::round((Double_t)nEntries/20.);

	TTreeReader _ntuple_reader(_ntuple_chain);
	TTreeReaderArrayValue<Float_t> _etaSCEle(_ntuple_reader, "etaSCEle");
	TTreeReaderArrayValue<Float_t> _phiSCEle(_ntuple_reader, "phiSCEle");
	TTreeReaderArrayValue<Float_t> _R9Ele(_ntuple_reader, "R9Ele");
	TTreeReaderArrayValue<Float_t> _pAtVtxGsfEle(_ntuple_reader, "pAtVtxGsfEle");
	TTreeReaderAnyValue<Float_t> _invMass_ECAL_ele(_ntuple_reader, "invMass_ECAL_ele");

	TH2D phi_vs_eta("phi_vs_eta_anom", "\\eta^{SC} \\ vs \\ \\phi^{SC}; \\phi^{SC}; \\eta^{SC}", 63, -3.15, 3.15, 120, -3.0, 3.0);
	TH1F R9dist("R9dist", "R_{9}; R_{9}; # of Z #rightarrow ee events", 300, -1.5, 1.5);

	while (_ntuple_reader.Next()) {
		ULong64_t current_entry = _ntuple_reader.GetCurrentEntry();
		if(current_entry % reportEvery == 0 ) {
			std::cout<<"\tProcessing entry "<<current_entry<<" ("<< (Double_t)current_entry*100./(Double_t)nEntries<<" %) "<<std::endl
			<<"\t\t"<< _ntuple_chain->GetFile()->GetName()<<std::endl;
		}

		Float_t invMass = _invMass_ECAL_ele;
		if(invMass < 70. || invMass > 110.) continue;

		Char_t passedPtCut = 0;

		for(auto iEle : INDICES){
			Float_t pT_ele = _pAtVtxGsfEle[iEle];
			if(pT_ele > 30.) setBit(passedPtCut, iEle, 1);
		}

		if(passedPtCut != 3) continue;

		for(auto iEle : INDICES){
			Float_t R9Ele = _R9Ele[iEle];
			R9dist.Fill(R9Ele);
			if( R9Ele < 1.) continue;
			Float_t etaSCEle = _etaSCEle[iEle];
			Float_t phiSCEle = _phiSCEle[iEle];
			phi_vs_eta.Fill(phiSCEle, etaSCEle);
		}
	}


	TFile *file = _ntuple_chain->GetCurrentFile();
	_ntuple_chain->SetDirectory(0);
	file->Delete();


	TFile outRoot(_out_file.c_str(), "RECREATE");
	outRoot.cd();
	phi_vs_eta.Write();
	R9dist.Write();
	outRoot.Close();
	std::cout<<"Graphs written to "<<_out_file<<std::endl;
};