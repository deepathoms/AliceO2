// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD Raw data format generation                                           //
//  Take the output of the trapsimulator the [2-4]x32 bit words and          //
//  associated headers and links etc. and produce the output of the cru.     //
//  Hence the incredibly original name.                                      //
//                                                                           //
//  Some parts bare a striking resemblence to (at least when written)        //
//  FIT/FTO/simulation/src/digit2raw.cxx                                     //
///////////////////////////////////////////////////////////////////////////////

#include "fairlogger/Logger.h"

#include <DataFormatsTRD/RawData.h>
#include <DataFormatsTRD/Tracklet64.h>
#include <DataFormatsTRD/LinkRecord.h>
#include <DataFormatsTRD/TriggerRecord.h>

#include "CommonUtils/StringUtils.h"
#include "CommonUtils/ConfigurableParam.h"
#include "DetectorsCommonDataFormats/NameConf.h"
#include "DetectorsRaw/HBFUtils.h"
#include "TRDSimulation/Trap2CRU.h"
#include "DataFormatsParameters/GRPObject.h"

#include <boost/program_options.hpp>
#include <TSystem.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <string>
#include <iomanip>
#include <iostream>
#include <iomanip>
#include "TCanvas.h"
#include <TTree.h>
#include <TFile.h>
#include <ostream>
#include <fstream>

namespace bpo = boost::program_options;

void trap2raw(const std::string& inpName, const std::string& outDir, int verbosity,
              bool filePerLink, uint32_t rdhV = 4, bool noEmptyHBF = false,
              int superPageSizeInB = 1024 * 1024);

int main(int argc, char** argv)
{
  bpo::variables_map vm;
  bpo::options_description opt_general("Usage:\n  " + std::string(argv[0]) +
                                       "Convert TRD sim otuput to CRU raw data\n");
  bpo::options_description opt_hidden("");
  bpo::options_description opt_all;
  bpo::positional_options_description opt_pos;

  try {
    auto add_option = opt_general.add_options();
    add_option("help,h", "Print this help message");
    add_option("verbosity,v", bpo::value<int>()->default_value(0), "verbosity level");
    add_option("input-file,i", bpo::value<std::string>()->default_value("trdtrapraw.root"), "input Trapsim raw file");
    add_option("file-per-halfcru,l", bpo::value<bool>()->default_value(false)->implicit_value(true), "create output file per half cru, 2 files per cru 15 links each.");
    add_option("output-dir,o", bpo::value<std::string>()->default_value("./"), "output directory for raw data");
    uint32_t defRDH = o2::raw::RDHUtils::getVersion<o2::header::RAWDataHeader>();
    add_option("rdh-version,r", bpo::value<uint32_t>()->default_value(defRDH), "RDH version to use");
    add_option("no-empty-hbf,e", bpo::value<bool>()->default_value(false)->implicit_value(true), "do not create empty HBF pages (except for HBF starting TF)");
    //    add_option("data-format,d", bpo::value<uint32_t>()->default_value(1)->implicit_value(true), "format of data, default 1, see documentation.");
    add_option("configKeyValues", bpo::value<std::string>()->default_value(""), "comma-separated configKeyValues");

    opt_all.add(opt_general).add(opt_hidden);
    bpo::store(bpo::command_line_parser(argc, argv).options(opt_all).positional(opt_pos).run(), vm);

    if (vm.count("help")) {
      std::cout << opt_general << std::endl;
      exit(0);
    }

    bpo::notify(vm);
  } catch (bpo::error& e) {
    std::cerr << "ERROR: " << e.what() << std::endl
              << std::endl;
    std::cerr << opt_general << std::endl;
    exit(1);
  } catch (std::exception& e) {
    std::cerr << e.what() << ", application will now exit" << std::endl;
    exit(2);
  }
  //  o2::conf::ConfigurableParam::updateFromString(vm["configKeyValues"].as<std::string>());

  std::cout << "yay it ran" << std::endl;
  trap2raw(vm["input-file"].as<std::string>(), vm["output-dir"].as<std::string>(), vm["verbosity"].as<int>(),
           vm["file-per-halfcru"].as<bool>(), vm["rdh-version"].as<uint32_t>(), vm["no-empty-hbf"].as<bool>());

  return 0;
}

void trap2raw(const std::string& inpName, const std::string& outDir, int verbosity, bool filePerLink, uint32_t rdhV, bool noEmptyHBF, int superPageSizeInB)
{

  TStopwatch swTot;
  swTot.Start();
  o2::trd::Trap2CRU mc2raw;
  //mc2raw.setFilePerLink(filePerLink); TODO come back to this, lets get it working first.
  mc2raw.setVerbosity(verbosity);
  auto& wr = mc2raw.getWriter();
  std::string inputGRP = o2::base::NameConf::getGRPFileName();
  const auto grp = o2::parameters::GRPObject::loadFrom(inputGRP);
  // wr.setContinuousReadout(grp->isDetContinuousReadOut(o2::detectors::DetID::TRD)); // must be set explicitly
  wr.setContinuousReadout(false); // above should work, but I know this is correct, come back TODO
  wr.setSuperPageSize(superPageSizeInB);
  wr.useRDHVersion(rdhV);
  wr.setDontFillEmptyHBF(noEmptyHBF);

  std::string outDirName(outDir);
  if (outDirName.back() != '/') {
    outDirName += '/';
  }
  // if needed, create output directory
  if (gSystem->AccessPathName(outDirName.c_str())) {
    if (gSystem->mkdir(outDirName.c_str(), kTRUE)) {
      LOG(FATAL) << "could not create output directory " << outDirName;
    } else {
      LOG(INFO) << "created output directory " << outDirName;
    }
  }

  mc2raw.readTrapData(outDirName, inpName, superPageSizeInB);
  wr.writeConfFile(wr.getOrigin().str, "RAWDATA", o2::utils::concat_string(outDirName, wr.getOrigin().str, "raw.cfg"));
  //
  swTot.Stop();
  swTot.Print();
}
