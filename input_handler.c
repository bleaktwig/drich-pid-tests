#include <cstdlib>
#include <stdio.h>

// C++.
#include <iostream>
#include <vector>
#include <set>

// ROOT.
#include "TSystem.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRegexp.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TVector3.h"
#include "TLeaf.h"
#include "TParameter.h"
#include "TObjString.h"
#include "TParticle.h"
#include <ROOT/RVec.hxx>

// EDM4HEP and EDM4EIC.
R__LOAD_LIBRARY(edm4hep)
R__LOAD_LIBRARY(edm4eic)
R__LOAD_LIBRARY(edm4hepDict)
R__LOAD_LIBRARY(edm4eicDict)

// DD4HEP.
#include <DD4hep/Detector.h>
#include <DDRec/DetectorData.h>
#include "DDRec/CellIDPositionConverter.h"

// local
// #include "WhichRICH.h"

TTree* eicreconTree = NULL;

TFile* outFile = NULL;
TTree* outTree = NULL;

float mc_p;
float mc_eta;
float mc_phi;

int   rec_track;
float track_p;
float track_eta;
float track_phi;

float h_endcap_clust_E;
float h_endcap_clust_eta;
float h_endcap_clust_phi;

float e_endcap_clust_E;
float e_endcap_clust_eta;
float e_endcap_clust_phi;

float barrel_clust_E;
float barrel_clust_eta;
float barrel_clust_phi;

void SetInputBranchAddresses();
void CreateOutputTree(TString outFileName);
void ResetVariables();

using namespace std;

/**
 * Constants to go from cell_id to local coordinate system. Formula is:
 *     local_x = DILATION * pos_x + OFFSET
 *     local_y = DILATION * pos_y + OFFSET
 */
#define DILATION 4.3
#define OFFSET   0.5

/**
 * Print an array of basic information for every pixel of each pdu in sector 0
 *     to a csv file. The data includes the pixel's parameters (sector, pdu,
 *     sipm, x, y). Then, print the cell's position in the local coordinate sys-
 *     tem (lx, ly).
 */
int get_cells_data(
    const dd4hep::DetElement dRICH,
    const dd4hep::DDSegmentation::BitFieldCoder *readoutCoder
) {
    // Open file.
    FILE *fout = fopen("csv/cells.csv", "w");

    // Print the header.
    fprintf(fout, "sector,pdu,sipm,x,y,lx,ly\n");

    // Iterate through the children of dRICH.
    for (auto const &[d_name, d_sensor] : dRICH.children()) {
        if (d_name.find("sensor_de_sec0") == std::string::npos) continue;
        const auto *det_pars =
            d_sensor.extension<dd4hep::rec::VariantParameters>(true);

        // Get sensor ID.
        ULong_t cell_id = ULong_t(d_sensor.id());

        // Get data from sensor ID.
        uint64_t sector = readoutCoder->get(cell_id, "sector"); //  3 bits.
        uint64_t pdu    = readoutCoder->get(cell_id, "pdu");    // 12 bits.
        uint64_t sipm   = readoutCoder->get(cell_id, "sipm");   //  6 bits.
        uint64_t x      = readoutCoder->get(cell_id, "x");      // 16 bits.
        uint64_t y      = readoutCoder->get(cell_id, "y");      // 16 bits.

        // Get position in local coordinate system.
        double pixel_x = DILATION*det_pars->get<double>("pos_x") + OFFSET;
        double pixel_y = DILATION*det_pars->get<double>("pos_y") + OFFSET;

        // Print the set of 64 pixels.
        for (int xi = 0; xi < 8; ++xi) {
            for (int yi = 0; yi < 8; ++yi) {
                fprintf(
                    fout,
                    "%lu,%lu,%lu,%lu,%lu,%.2f,%.2f\n",
                    sector, pdu, sipm, x+xi, y+yi, pixel_x+xi, pixel_y+yi
                );
            }
        }
    }

    fclose(fout);

    return 0;
}

int extractSimuReco(TString inFileSimu, TString inFileReco) {
    // Setup detector instance.
    // NOTE. This assumes that the user is running on the eic-shell, with the
    //       dRICH environ.sh sourced.
    dd4hep::Detector *det = &(dd4hep::Detector::getInstance());
    det->fromXML("/opt/detector/epic-23.10.0/share/epic/epic.xml");
    const dd4hep::DetElement dRICH = det->detector("DRICH");

    // Get BitFieldCoder to decode the cellID.
    const dd4hep::DDSegmentation::BitFieldCoder *readoutCoder =
        det->readout("DRICHHits").idSpec().decoder();

    // Get CellIDPositionConverter to get the pixel's position in the ePIC glo-
    //     bal coordinate system.
    dd4hep::rec::CellIDPositionConverter geoConverter(*det);

    // Get TTreeRReader for reconstructed file.
    TFile *recoFile = new TFile(inFileReco);
    TTree *recoTree = (TTree *) recoFile->Get("events");
    TTreeReader recoTR(recoTree);

    // Get cell_id.
    TTreeReaderArray<uint64_t> RH_cellID(recoTR, "DRICHRawHits.cell_id");
    TTreeReaderArray<int32_t>  RH_charge(recoTR, "DRICHRawHits.charge");
    TTreeReaderArray<int32_t>  RH_time(  recoTR, "DRICHRawHits.timeStamp");

    // Set treereader to first entry.
    recoTR.SetEntry(-1);

    // -------------------------------------------------------------------------
    get_cells_data(dRICH, readoutCoder);

    // -------------------------------------------------------------------------
    // Print all hits in .csv format.
    // printf("event,sector,pdu,sipm,lx,ly,x,y,z,charge,time\n");
    // uint64_t event_i = -1;
    // while(recoTR.Next()) {
    //     ++event_i;
    //     int nhits = RH_cellID.GetSize();
    //     for (int hit_i = 0; hit_i < nhits; ++hit_i) {
    //         uint64_t cell_id = RH_cellID[hit_i];
    //
    //         uint64_t sector = readoutCoder->get(cell_id, "sector"); //  3 bits.
    //         uint64_t pdu    = readoutCoder->get(cell_id, "pdu");    // 12 bits.
    //         uint64_t sipm   = readoutCoder->get(cell_id, "sipm");   //  6 bits.
    //         uint64_t x      = readoutCoder->get(cell_id, "x");      // 16 bits.
    //         uint64_t y      = readoutCoder->get(cell_id, "y");      // 16 bits.
    //
    //         dd4hep::Position point = geoConverter.position(cell_id);
    //         printf(
    //             "%lu,%lu,%lu,%lu,%lu,%lu,%.2f,%.2f,%.2f,%u,%u\n",
    //             event_i, sector, pdu, sipm, x, y, point.x(), point.y(), point.z(),
    //             RH_charge[hit_i], RH_time[hit_i]
    //         );
    //     }
    // }
    // -------------------------------------------------------------------------
    // Print all cell numbers and positions in .csv format.
    // printf("sector,pdu,sipm,lx,ly,x,y,z\n");
    // for (uint64_t sec = 0; sec < 6; ++sec) {
    //     for (uint64_t pdu = 0; pdu < 277; ++pdu) {
    //         for (uint64_t sipm = 0; sipm < 4; ++sipm) {
    //             for (uint64_t x = 0; x < 8; ++x) {
    //                 for (uint64_t y = 0; y < 8; ++y) {
    //                     uint64_t cell_id = make_cellID(
    //                         0b01111000, sec, pdu, sipm, x, y
    //                     );
    //                     dd4hep::Position point = geoConverter.position(cell_id);
    //                     printf(
    //                         "%lu,%lu,%lu,%lu,%lu,%.2f,%.2f,%.2f\n",
    //                         sec, pdu, sipm, x, y, point.x(), point.y(), point.z()
    //                     );
    //
    //                 }
    //             }
    //         }
    //     }
    // }
    // -------------------------------------------------------------------------

    return 0;
}

int input_handler() {
    extractSimuReco(
        "/home/twig/code/eic/drich-dev/out/sim.edm4hep.root", // simu.
        "/home/twig/code/eic/drich-dev/out/rec.noise.edm4hep.root" // reco.
    );

    return 0;
}
