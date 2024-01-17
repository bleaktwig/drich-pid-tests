// C.
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// C++.
#include <map>

// ROOT.
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"

// dd4hep.
#include <DD4hep/Detector.h>
#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>

/**
 * Strings.
 */
#define CELLMAP_LOC "csv/cellmap.csv"
#define HITMAP_LOC  "csv/hitmap.csv"
#define BUFFERSIZE 256

/**
 * Constants to go from cell ID to local coordinate system. Formula is:
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
    FILE *f_map = fopen(CELLMAP_LOC, "w");

    // Write header.
    fprintf(f_map, "cell_id,x,y\n");

    // Iterate through the dRICH sub-detectors.
    for (auto const &[d_name, d_sensor] : dRICH.children()) {
        // NOTE. Since the coordinate system is the same for all sectors, we on-
        //       ly care about sector 0 here.
        if (d_name.find("sensor_de_sec0") == std::string::npos) continue;
        const auto *det_pars =
            d_sensor.extension<dd4hep::rec::VariantParameters>(true);

        // Get sensor ID.
        uint64_t cell_id = (uint64_t) d_sensor.id();

        // Set system to dRICH.
        readout_coder->set(cell_id, "system", 120); // TODO. Remove hardcoding.

        // Get position in local coordinate system.
        double pixel_x = DILATION*det_pars->get<double>("pos_x") + OFFSET;
        double pixel_y = DILATION*det_pars->get<double>("pos_y") + OFFSET;

        // Get matrix index from pixel position.
        uint64_t idx_x = ((uint64_t) pixel_x) - PIXELX_MIN + 1;
        uint64_t idx_y = ((uint64_t) pixel_y) - PIXELY_MIN + 1;

        // Walk through the set of 64 pixels.
        for (int xi = 0; xi < 8; ++xi) {
            for (int yi = 0; yi < 8; ++yi) {
                // Include xi and yi into the cell ID.
                readout_coder->set(cell_id, "x", xi);
                readout_coder->set(cell_id, "y", yi);

                // Print into the output file.
                fprintf(f_map, "%lu,%lu,%lu\n", cell_id, idx_x+xi, idx_y+yi);
            }
        }
    }

    // Clean-up.
    fclose(f_map);

    return 0;
}

/**
 * Read cellmap from CELLMAP_LOC into an std::map. If CELLMAP_LOC doesn't exist,
 *     create it and then read it.
 */
int read_cellmap(
    const dd4hep::DetElement dRICH,
    const dd4hep::DDSegmentation::BitFieldCoder *readout_coder,
    std::map<uint64_t, std::pair<uint64_t, uint64_t>> *cellmap
) {
    // Open (or create and open) cellmap.
    FILE *f_map = fopen(CELLMAP_LOC, "r");
    if (!f_map) {
        create_cellmap(dRICH, readout_coder);
        f_map = fopen(CELLMAP_LOC, "r");
    }

    // Iterate through cellmap csv.
    char buffer[BUFFERSIZE];
    fgets(buffer, BUFFERSIZE, f_map); // Ignore header.
    while (fgets(buffer, BUFFERSIZE, f_map)) {

        // Rescue each number from the cellmap file.
        uint64_t data[3] = {0, 0, 0};
        char *val = strtok(buffer, ",");
        for (int val_i = 0; val != NULL; ++val_i) {
            data[val_i] = strtoull(val, NULL, 0);
            val = strtok(NULL, ",");
        }

        // Write data into map.
        (*cellmap)[data[0]] = {data[1], data[2]};
    }

    // Clean-up.
    fclose(f_map);

    return 0;
}

/**
 * Extract hits from a simulated and reconstructed ROOT file, associate them,
 *     and write them to a csv file.
 */
int extractSimuReco(TString f_sim, TString f_rec) {
    // Setup detector instance.
    // NOTE. This assumes that the user is running on the eic-shell, with the
    //       dRICH environ.sh sourced.
    dd4hep::Detector *det = &(dd4hep::Detector::getInstance());
    det->fromXML("/opt/detector/epic-23.10.0/share/epic/epic.xml");
    const dd4hep::DetElement dRICH = det->detector("DRICH");

    // Get BitFieldCoder to decode the cellID.
    const dd4hep::DDSegmentation::BitFieldCoder *readout_coder =
        det->readout("DRICHHits").idSpec().decoder();

    // Read a cell map from CELLMAP_LOC.
    std::map<uint64_t, std::pair<uint64_t, uint64_t>> cellmap;
    read_cellmap(dRICH, readout_coder, &cellmap);

    // Get sim and rec TTreeReaders.
    TTreeReader s_tree((TTree *) (new TFile(f_sim))->Get("events"));
    TTreeReader r_tree((TTree *) (new TFile(f_rec))->Get("events"));

    // Associate TTreeReaderArrays with relevant data from trees.
    // Sim.
    TTreeReaderArray<uint64_t> s_cell_id(s_tree, "DRICHHits.cellID");
    TTreeReaderArray<int>      s_index(  s_tree, "DRICHHits#0.index");

    // Rec.
    TTreeReaderArray<uint64_t> r_cell_id(r_tree, "DRICHRawHits.cellID");
    TTreeReaderArray<int32_t>  r_charge( r_tree, "DRICHRawHits.charge");
    TTreeReaderArray<int32_t>  r_time(   r_tree, "DRICHRawHits.timeStamp");

    // Set TTreeReaders to first entry.
    s_tree.SetEntry(-1);
    r_tree.SetEntry(-1);

    // Open hitmap output file.
    FILE *f_map = fopen(HITMAP_LOC, "w");

    // Print header.
    // printf("sector,pdu,sipm,x,y,time,charge,pindex\n");
    fprintf(f_map, "event,time,sector,x,y,charge,pindex\n");

    // Iterate through events.
    uint64_t event_i = 0;
    while(s_tree.Next() && r_tree.Next()) {
        // NOTE. It might be a good idea to put hits into a list, order them by
        //       time, and then store them in the output file.

        // Iterate through rec and sim hits.
        for (int rh_it = 0; rh_it < r_cell_id.GetSize(); ++rh_it) {
            for (int sh_it = 0; sh_it < s_cell_id.GetSize(); ++sh_it) {
                if (s_cell_id[sh_it] != r_cell_id[rh_it]) continue;
                uint64_t cell_id = r_cell_id[rh_it];

                // Store sector and check hit position from cellmap.
                uint64_t sector = readout_coder->get(cell_id, "sector");
                readout_coder->set(cell_id, "sector", 0);

                if (!cellmap.count(cell_id)) {
                    printf("Cell ID %lu not in cellmap. Exiting.\n", cell_id);
                    return 1;
                }

                // Write to stdout.
                fprintf(
                    f_map, "%lu,%d,%lu,%lu,%lu,%d,%d\n",
                    event_i, r_time[rh_it], sector,
                    cellmap[cell_id].first, cellmap[cell_id].second,
                    r_charge[rh_it], s_index[sh_it]
                );
            }
        }
        ++event_i;
    }

    // Clean-up.
    fclose(f_map);

    return 0;
}

int input_handler() {
    extractSimuReco(
        "/home/twig/code/eic/drich-dev/out/sim.edm4hep.root", // simu.
        "/home/twig/code/eic/drich-dev/out/rec.noise.edm4hep.root" // reco.
    );

    return 0;
}
