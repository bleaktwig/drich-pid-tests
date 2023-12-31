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
    const dd4hep::DDSegmentation::BitFieldCoder *readout_coder
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
        uint64_t sector = readout_coder->get(cell_id, "sector"); //  3 bits.
        uint64_t pdu    = readout_coder->get(cell_id, "pdu");    // 12 bits.
        uint64_t sipm   = readout_coder->get(cell_id, "sipm");   //  6 bits.
        uint64_t x      = readout_coder->get(cell_id, "x");      // 16 bits.
        uint64_t y      = readout_coder->get(cell_id, "y");      // 16 bits.

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

    // Get CellIDPositionConverter to get the pixel's position in the ePIC glo-
    //     bal coordinate system.
    // dd4hep::rec::CellIDPositionConverter geo_converter(*det);

    // Get TTreeReaders from simulated and reconstructed file.
    TTreeReader s_tree((TTree *) (new TFile(fsimu))->Get("events"));
    TTreeReader r_tree((TTree *) (new TFile(freco))->Get("events"));

    // Associate TTreeReaderArrays with relevant data from trees.
    TTreeReaderArray<uint64_t> r_cell_id(r_tree, "DRICHRawHits.cellID");
    TTreeReaderArray<int32_t>  r_charge( r_tree, "DRICHRawHits.charge");
    TTreeReaderArray<int32_t>  r_time(   r_tree, "DRICHRawHits.timeStamp");

    // Set TTreeReaders to first entry.
    s_tree.SetEntry(-1);
    r_tree.SetEntry(-1);

    // -------------------------------------------------------------------------
    // Print header.
    printf("sector,x,y,time,charge\n");

    // Iterate through events.
    while(r_tree.Next()) {
        // Iterate through hits.
        int nhits = r_cell_id.GetSize();
        for (int hit_i = 0; hit_i < nhits; ++hit_i) {
            uint64_t sector = readout_coder->get(r_cell_id[hit_i], "sector");
            uint64_t pdu    = readout_coder->get(r_cell_id[hit_i], "pdu");
            uint64_t sipm   = readout_coder->get(r_cell_id[hit_i], "sipm");
            uint64_t x      = readout_coder->get(r_cell_id[hit_i], "x");
            uint64_t y      = readout_coder->get(r_cell_id[hit_i], "y");

            // TODO. Find a way to associate cellID (or its components) to a
            //       position in the local coordinate system. Then print that.
            //       That's all folks!
        }
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
