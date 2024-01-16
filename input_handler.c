// C.
#include <cstdlib>
#include <stdio.h>

// ROOT.
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"

// dd4hep.
#include <DD4hep/Detector.h>
#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>

/**
 * Constants to go from cell_id to local coordinate system. Formula is:
 *     local_x = DILATION * pos_x + OFFSET
 *     local_y = DILATION * pos_y + OFFSET
 */
#define DILATION 4.3
#define OFFSET   0.5

/**
 * Position boundaries in the local coordinate system. Based on these numbers,
 * the input matrix size is 277x609, and thus has 168693 pixels.
 *
 * TODO. Store these when generating cells.csv instead of setting them as
 *       constants.
 */
#define PIXELX_MIN  485
#define PIXELX_MAX  761
#define PIXELY_MIN -304
#define PIXELY_MAX  304

/**
 * Print a map from cell ID to matrix index to a csv file. This function assumes
 *     that DILATION, OFFSET, PIXELX_MIN, and PIXELY_MIN are correctly setup to
 *     the detector setup.
 */
int create_cellmap(
    const dd4hep::DetElement dRICH,
    const dd4hep::DDSegmentation::BitFieldCoder *readout_coder
) {
    // Create cellmap file.
    FILE *fout = fopen("csv/cellmap.csv", "w");

    // Iterate through the dRICH sub-detectors.
    for (auto const &[d_name, d_sensor] : dRICH.children()) {
        // NOTE. Since the coordinate system is the same for all sectors, we on-
        //       ly care about sector 0 here.
        if (d_name.find("sensor_de_sec0") == std::string::npos) continue;
        const auto *det_pars =
            d_sensor.extension<dd4hep::rec::VariantParameters>(true);

        // Get sensor ID.
        uint64_t cell_id = (uint64_t) d_sensor.id();

        // Decode relevant variables in sensor ID.
        uint64_t system = readout_coder->get(cell_id, "system"); //  8 bits.
        uint64_t sector = readout_coder->get(cell_id, "sector"); //  3 bits.
        uint64_t pdu    = readout_coder->get(cell_id, "pdu");    // 12 bits.
        uint64_t sipm   = readout_coder->get(cell_id, "sipm");   //  6 bits.

        // Get position in local coordinate system.
        double pixel_x = DILATION*det_pars->get<double>("pos_x") + OFFSET;
        double pixel_y = DILATION*det_pars->get<double>("pos_y") + OFFSET;

        // Get matrix index from pixel position.
        uint64_t idx_x = ((uint64_t) pixel_x) - PIXELX_MIN + 1;
        uint64_t idx_y = ((uint64_t) pixel_y) - PIXELY_MIN + 1;

        // Print the set of 64 pixels to the output file.
        for (int xi = 0; xi < 8; ++xi) {
            for (int yi = 0; yi < 8; ++yi) {
                // Include xi and yi into the cell ID.
                readout_coder->set(cell_id, "x", xi);
                readout_coder->set(cell_id, "y", yi);
                fprintf(fout, "%lu,%lu,%lu\n", cell_id, idx_x+xi, idx_y+yi);
            }
        }
    }

    fclose(fout);

    return 0;
}

/**
 * Extract hits from a simulated and reconstructed ROOT file, associate them,
 *     and write them to a csv file.
 */
int extractSimuReco(TString fsimu, TString freco) {
    // Setup detector instance.
    // NOTE. This assumes that the user is running on the eic-shell, with the
    //       dRICH environ.sh sourced.
    dd4hep::Detector *det = &(dd4hep::Detector::getInstance());
    det->fromXML("/opt/detector/epic-23.10.0/share/epic/epic.xml");
    const dd4hep::DetElement dRICH = det->detector("DRICH");

    // Get BitFieldCoder to decode the cellID.
    const dd4hep::DDSegmentation::BitFieldCoder *readout_coder =
        det->readout("DRICHHits").idSpec().decoder();

    create_cellmap(dRICH, readout_coder);
    return 0;

    // Get CellIDPositionConverter to get the pixel's position in the ePIC glo-
    //     bal coordinate system.
    // dd4hep::rec::CellIDPositionConverter geo_converter(*det);

    // Get TTreeReaders from simulated and reconstructed file.
    TTreeReader s_tree((TTree *) (new TFile(fsimu))->Get("events"));
    TTreeReader r_tree((TTree *) (new TFile(freco))->Get("events"));

    // Associate TTreeReaderArrays with relevant data from trees.
    // Simu.
    TTreeReaderArray<uint64_t> s_cell_id(s_tree, "DRICHHits.cellID");
    TTreeReaderArray<int>      s_index(  s_tree, "DRICHHits#0.index");

    // Reco.
    TTreeReaderArray<uint64_t> r_cell_id(r_tree, "DRICHRawHits.cellID");
    TTreeReaderArray<int32_t>  r_charge( r_tree, "DRICHRawHits.charge");
    TTreeReaderArray<int32_t>  r_time(   r_tree, "DRICHRawHits.timeStamp");

    // Set TTreeReaders to first entry.
    s_tree.SetEntry(-1);
    r_tree.SetEntry(-1);

    // -------------------------------------------------------------------------
    // Print header.
    printf("sector,pdu,sipm,x,y,time,charge,pindex\n");

    // Iterate through events.
    while(s_tree.Next() && r_tree.Next()) {
        // Iterate through reconstructed hits.
        for (int rhit_i = 0; rhit_i < r_cell_id.GetSize(); ++rhit_i) {
            // Iterate through simulated hits.
            for (int shit_i = 0; shit_i < s_cell_id.GetSize(); ++shit_i) {
                // Only continue if we're working with the same hit.
                if (s_cell_id[shit_i] != r_cell_id[rhit_i]) continue;

                // Decode hit from cellID.
                uint64_t sector = readout_coder->get(r_cell_id[rhit_i], "sector");
                uint64_t pdu    = readout_coder->get(r_cell_id[rhit_i], "pdu");
                uint64_t sipm   = readout_coder->get(r_cell_id[rhit_i], "sipm");
                uint64_t x      = readout_coder->get(r_cell_id[rhit_i], "x");
                uint64_t y      = readout_coder->get(r_cell_id[rhit_i], "y");

                printf(
                    "%lu,%lu,%lu,%lu,%lu,%d,%d,%d\n",
                    sector, pdu, sipm, x, y,
                    r_time[rhit_i], r_charge[rhit_i],
                    s_index[shit_i]
                );


            }
        }

        break;
    }

    return 0;
}

int input_handler() {
    extractSimuReco(
        "/home/twig/code/eic/drich-dev/out/sim.edm4hep.root", // simu.
        "/home/twig/code/eic/drich-dev/out/rec.noise.edm4hep.root" // reco.
    );

    return 0;
}
