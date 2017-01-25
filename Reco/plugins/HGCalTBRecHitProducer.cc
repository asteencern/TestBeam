#include "HGCal/Reco/plugins/HGCalTBRecHitProducer.h"
#include <iostream>

HGCalTBRecHitProducer::HGCalTBRecHitProducer(const edm::ParameterSet& cfg)
	: outputCollectionName(cfg.getParameter<std::string>("OutputCollectionName")),
	  _digisToken(consumes<HGCalTBDigiCollection>(cfg.getParameter<edm::InputTag>("digiCollection"))),
	  _pedestalLow_filename(cfg.getParameter<std::string>("pedestalLow")),
	  _pedestalHigh_filename(cfg.getParameter<std::string>("pedestalHigh")),
	  _gainsLow_filename(cfg.getParameter<std::string>("gainLow")),
	  _gainsHigh_filename(cfg.getParameter<std::string>("gainHigh")),
	  _adcSaturation(cfg.getParameter<int>("adcSaturation")),
	  //_LG2HG_value(cfg.getParameter<std::vector<double> >("LG2HG_CERN")),
	  //_mapFile(cfg.getParameter<std::string>("mapFile")),
	  _layers_config(cfg.getParameter<int>("layers_config")),
	  _commonModeThreshold(cfg.getUntrackedParameter<double>("CommonModeThreshold",2)), // in MIP
	  _doCommonMode(cfg.getUntrackedParameter<bool>("doCommonMode",true)),
	  _mipToGeV(cfg.getUntrackedParameter<double>("MIPToMeV",52.81e-06)) //obtained from simulation
{
	produces <HGCalTBRecHitCollection>(outputCollectionName);
	//	std::cout << " >>> _LG2HG_value size = " << _LG2HG_value.size() << std::endl;

	if(_layers_config == 0){
	  _LG2HG_value = cfg.getParameter<std::vector<double> >("LG2HG_FNAL");
	  _mapFile = cfg.getParameter<std::string>("mapFile_FNAL");
	}
	else{
	  _LG2HG_value = cfg.getParameter<std::vector<double> >("LG2HG_CERN");
	  _mapFile = cfg.getParameter<std::string>("mapFile_CERN");
	}
	
	std::cout << cfg.dump() << std::endl;

	HGCalCondObjectTextIO io(0);
	edm::FileInPath fip(_mapFile);
        if (!io.load(fip.fullPath(), essource_.emap_)) {
	  throw cms::Exception("Unable to load electronics map");
	};
	
	if( _doCommonMode )
	  rhcm = new RecHitCommonMode( essource_.emap_ );

	double adctomip[]={15.926 ,15.6252 ,14.8598 ,15.9285 ,16.873 ,15.6593 ,15.6549 ,14.5163 ,16.2659 ,15.5812 ,15.369 ,16.4178 ,14.4392 ,15.2915 ,15.51 ,14.2051};
	std::vector<double> vec; vec.insert( vec.begin(), adctomip, adctomip+16 );
	_skirocADCToMip = cfg.getUntrackedParameter< std::vector<double> >("skirocADCToMip",vec);
	_commonModeThreshold=_commonModeThreshold*_mipToGeV;
}

HGCalTBRecHitProducer::~HGCalTBRecHitProducer()
{
  if( _doCommonMode )
    delete rhcm;
}

void HGCalTBRecHitProducer::produce(edm::Event& event, const edm::EventSetup& iSetup)
{

	std::auto_ptr<HGCalTBRecHitCollection> rechits(new HGCalTBRecHitCollection); //auto_ptr are automatically deleted when out of scope

	edm::Handle<HGCalTBDigiCollection> digisHandle;
	event.getByToken(_digisToken, digisHandle);

#ifdef ESPRODUCER
// the conditions should be defined by ESProduers and available in the iSetup
	edm::ESHandle<HGCalCondPedestals> pedestalsHandle;
	//iSetup.get<HGCalCondPedestals>().get(pedestalsHandle); ///\todo should we support such a syntax?
	HGCalCondPedestals pedestals = iSetup.get<HGCalCondPedestals>();

	edm::ESHandle<HGCalCondGains>     adcToGeVHandle;
	//iSetup.get<HGCalCondGains>().get(adcToGeVHandle);
	HGCalCondGains adcToGeV = iSetup.get<HGCalCondGains>();
#else
	HGCalCondGains adcToGeV_low, adcToGeV_high;
	HGCalCondPedestals pedestals_low, pedestals_high;
	HGCalCondObjectTextIO condIO(HGCalTBNumberingScheme::scheme());
	//HGCalElectronicsMap emap;
	//    assert(io.load("mapfile.txt",emap)); ///\todo to be trasformed into exception
	if(_pedestalLow_filename != "") assert(condIO.load(_pedestalLow_filename, pedestals_low));
	if(_pedestalHigh_filename != "") 	assert(condIO.load(_pedestalHigh_filename, pedestals_high));
	if(_gainsLow_filename != "") 	assert(condIO.load(_gainsLow_filename, adcToGeV_low));
	if(_gainsHigh_filename != "") 	assert(condIO.load(_gainsHigh_filename, adcToGeV_high));
	///\todo check if reading the conditions from file some channels are not in the file!
#endif
	HGCalTBRecHitCollection tmp;
	for(auto digi_itr = digisHandle->begin(); digi_itr != digisHandle->end(); ++digi_itr) {
#ifdef DEBUG
		std::cout << "[RECHIT PRODUCER: digi]" << *digi_itr << std::endl;
#endif

//------------------------------ this part should go into HGCalTBRecoAlgo class

		SKIROC2DataFrame digi(*digi_itr);
		unsigned int nSamples = digi.samples();

		// if there are more than 1 sample, we need to define a reconstruction algorithm
		// now taking the first sample
		for(unsigned int iSample = 0; iSample < nSamples; ++iSample) {

			float pedestal_low_value = (pedestals_low.size() > 0) ? pedestals_low.get(digi.detid())->value : 0;
			float pedestal_high_value = (pedestals_high.size() > 0) ? pedestals_high.get(digi.detid())->value : 0;
			float adcToGeV_low_value = (adcToGeV_low.size() > 0) ? adcToGeV_low.get(digi.detid())->value : 1;
			float adcToGeV_high_value = (adcToGeV_high.size() > 0) ? adcToGeV_high.get(digi.detid())->value : 1;

			float energyLow = digi[iSample].adcLow() - pedestal_low_value * adcToGeV_low_value;
			float energyHigh = digi[iSample].adcHigh() - pedestal_high_value * adcToGeV_high_value;

			float energy = -1.;

			HGCalTBRecHit recHit(digi.detid(), energy, energyLow, energyHigh, digi[iSample].tdc()); //, _LG2HG_value, _gainThr_value * adcToGeV_high_value); ///\todo use time calibration!

			uint32_t EID = essource_.emap_.detId2eid(recHit.id());
			HGCalTBElectronicsId eid(EID);


			energy = ( energyHigh < _adcSaturation ) ? energyHigh : energyLow * _LG2HG_value.at(eid.iskiroc() - 1);
			recHit.setEnergy(energy/_skirocADCToMip[ eid.iskiroc()-1 ]*_mipToGeV);

			if(digi[iSample].adcHigh() > _adcSaturation) recHit.setFlag(HGCalTBRecHit::kHighGainSaturated);
			if(digi[iSample].adcLow() > _adcSaturation) recHit.setFlag(HGCalTBRecHit::kLowGainSaturated);

#ifdef DEBUG
			std::cout << recHit << std::endl;
#endif
			if(iSample == 0) tmp.push_back(recHit); ///\todo define an algorithm for the energy if more than 1 sample, code inefficient
		}
	}
	if( _doCommonMode ){
	  rhcm->evaluate( tmp,_commonModeThreshold );
	  for( auto hit:tmp ){
	    HGCalTBRecHit recHit(hit.id(), hit.energy()-rhcm->getMeanCommonModeNoise(hit.id()), hit.energyLow(), hit.energyHigh(), hit.time());
	    rechits->push_back(recHit);
	  }
	}
	else
	  for( auto hit:tmp ) rechits->push_back(hit);
	event.put(rechits, outputCollectionName);
}
// Should there be a destructor ??
DEFINE_FWK_MODULE(HGCalTBRecHitProducer);
