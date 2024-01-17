// C.
#include <limits.h>
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

/** dRICH system number for the cell ID. */
#define DRICH_SYSTEM 120

/** String handling constants. */
#define CELLMAP_LOC "csv/cellmap.csv"
#define BUFFERSIZE 256

/** Variables to go from global to local coordinates. */
#define DILATION 4.3
#define OFFSET   0.5

// Global instance of dRICH detector and bit field coder.
dd4hep::DetElement g_dRICH;
dd4hep::DDSegmentation::BitFieldCoder *g_readout_coder;

/**
 * Go from global coordinates to the dRICH sector local coordinate system.
 *
 * The local coordinate system is defined in such a way that each pixel becomes
 *     a flat square, flattening the entire dRICH sector. In this coordinate
 *     system, the dRICH rings become circles. The equation to go from a sen-
 *     sor's position in global coordinates to this system is, simply:
 *
 *         x_local = DILATION * x_global + OFFSET
 *         y_local = DILATION * y_global + OFFSET
 */
double global_to_local(double d) {return DILATION*d + OFFSET;}

/**
 * Initialize instances of dRICH detector and bit field coder from dd4hep. This
 *     function assumes that we are running over the eic-shell with the dRICH
 *     environ.sh script loaded.
 */
int init_dRICH() {
    dd4hep::Detector *det = &(dd4hep::Detector::getInstance());
    det->fromXML("/opt/detector/epic-23.10.0/share/epic/epic.xml");
    g_dRICH = det->detector("DRICH");
    g_readout_coder = det->readout("DRICHHits").idSpec().decoder();
    return 0;
}

/**
 * Print a cell map to a csv file under the name CELLMAP_LOC. The cell map index
 *     is the cell ID of each photomultiplier pixel in the sector 0 of dRICH.
 *     The stored data is the position of the pixel in the local coordinate sys-
 *     tem of each dRICH sector.
 */
int create_cellmap() {
    // Create cellmap file.
    FILE *f_map = fopen(CELLMAP_LOC, "w");

    // To define the 0,0 point of the matrix, find the lowest pixel position for
    //     x and y.
    int64_t pixel_x_min = INT_MAX;
    int64_t pixel_y_min = INT_MAX;

    for (auto const &[d_name, d_sensor] : g_dRICH.children()) {
        if (d_name.find("sensor_de_sec0") == std::string::npos) continue;
        const auto *det_pars =
            d_sensor.extension<dd4hep::rec::VariantParameters>(true);

        double pixel_x = global_to_local(det_pars->get<double>("pos_x"));
        double pixel_y = global_to_local(det_pars->get<double>("pos_y"));

        if (pixel_x < pixel_x_min) pixel_x_min = (int64_t) pixel_x;
        if (pixel_y < pixel_y_min) pixel_y_min = (int64_t) pixel_y;
    }

    // Write header.
    fprintf(f_map, "cell_id,x,y\n");

    // Iterate through the dRICH sub-detectors.
    for (auto const &[d_name, d_sensor] : g_dRICH.children()) {
        // NOTE. Since the coordinate system is the same for all sectors, we on-
        //       ly care about sector 0 here.
        if (d_name.find("sensor_de_sec0") == std::string::npos) continue;
        const auto *det_pars =
            d_sensor.extension<dd4hep::rec::VariantParameters>(true);

        // Get sensor ID.
        uint64_t cell_id = (uint64_t) d_sensor.id();

        // Set system to dRICH.
        g_readout_coder->set(cell_id, "system", DRICH_SYSTEM);

        // Get position in local coordinate system.
        double pixel_x = global_to_local(det_pars->get<double>("pos_x"));
        double pixel_y = global_to_local(det_pars->get<double>("pos_y"));

        // Get matrix index from pixel position.
        uint64_t idx_x = (uint64_t) (((int64_t) pixel_x) - pixel_x_min);
        uint64_t idx_y = (uint64_t) (((int64_t) pixel_y) - pixel_y_min);

        // Walk through the set of 64 pixels.
        for (int xi = 0; xi < 8; ++xi) {
            for (int yi = 0; yi < 8; ++yi) {
                // Include xi and yi into the cell ID.
                g_readout_coder->set(cell_id, "x", xi);
                g_readout_coder->set(cell_id, "y", yi);

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
int read_cellmap(std::map<uint64_t, std::pair<uint64_t, uint64_t>> *cellmap) {
    // Open (or create then open) cellmap.
    FILE *f_map = fopen(CELLMAP_LOC, "r");
    if (!f_map) {
        create_cellmap();
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
 * Write reconstructed hits into a csv file. Also writes simulated particle in-
 *     dex to each row for a many-to-one association.
 */
int write_rec_hits(const char *in_rec, const char *in_sim, const char *out) {
    // Read cell map from CELLMAP_LOC into an std::map.
    std::map<uint64_t, std::pair<uint64_t, uint64_t>> cellmap;
    read_cellmap(&cellmap);

    // Open trees.
    TTreeReader rec_tree((TTree *) (new TFile(in_rec))->Get("events"));
    TTreeReader sim_tree((TTree *) (new TFile(in_sim))->Get("events"));

    // Associate TTreeReaderArrays with relevant data from trees.
    TTreeReaderArray<uint64_t> r_cell_id(rec_tree, "DRICHRawHits.cellID");
    TTreeReaderArray<int32_t>  r_charge( rec_tree, "DRICHRawHits.charge");
    TTreeReaderArray<int32_t>  r_time(   rec_tree, "DRICHRawHits.timeStamp");

    TTreeReaderArray<uint64_t> s_cell_id(sim_tree, "DRICHHits.cellID");
    TTreeReaderArray<int>      s_index(  sim_tree, "DRICHHits#0.index");

    // Set TTreeReaders to first entry.
    sim_tree.SetEntry(-1);
    rec_tree.SetEntry(-1);

    // Open hitmap output file and print header.
    FILE *f_map = fopen(out, "w");
    fprintf(f_map, "event,time,sector,x,y,charge,pindex\n");

    // Iterate through events.
    uint64_t event_i = 0;
    while(sim_tree.Next() && rec_tree.Next()) {
        // NOTE. It might be a good idea to put hits into a list, order them by
        //       time, and then store them in the output file.

        // Iterate through rec and sim hits.
        for (int rh_it = 0; rh_it < r_cell_id.GetSize(); ++rh_it) {
            for (int sh_it = 0; sh_it < s_cell_id.GetSize(); ++sh_it) {
                if (s_cell_id[sh_it] != r_cell_id[rh_it]) continue;
                uint64_t cell_id = r_cell_id[rh_it];

                // Store sector and check hit position from cellmap.
                uint64_t sector = g_readout_coder->get(cell_id, "sector");
                g_readout_coder->set(cell_id, "sector", 0);

                if (!cellmap.count(cell_id)) {
                    fprintf(stderr, "%lu not in cellmap. Exiting.\n", cell_id);
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

/**
 * Write simulated particles into a csv file.
 */
int write_sim_parts(const char *in_sim, const char *out) {
    // NOTE. This could process only pindexes found in write_rec_hits, so as to
    //       avoid storing unnecessary information.

    // Open sim tree.
    TTreeReader sim_tree((TTree *) (new TFile(in_sim))->Get("events"));

    // Associate TTreeReaderArrays with relevant data from tree.
    TTreeReaderArray<int>    idx (sim_tree, "DRICHHits#0.index");
    TTreeReaderArray<int>    pdg (sim_tree, "MCParticles.PDG");
    TTreeReaderArray<float>  time(sim_tree, "MCParticles.time");
    TTreeReaderArray<double> vx  (sim_tree, "MCParticles.vertex.x");
    TTreeReaderArray<double> vy  (sim_tree, "MCParticles.vertex.y");
    TTreeReaderArray<double> vz  (sim_tree, "MCParticles.vertex.z");
    TTreeReaderArray<float>  px  (sim_tree, "MCParticles.momentum.x");
    TTreeReaderArray<float>  py  (sim_tree, "MCParticles.momentum.y");
    TTreeReaderArray<float>  pz  (sim_tree, "MCParticles.momentum.z");

    // Rewind TTreeReader.
    sim_tree.SetEntry(-1);

    // Open output file and print header.
    FILE *f_parts = fopen(out, "w");
    fprintf(f_parts, "event,pindex,pdg,time,vx,vy,vz,px,py,pz\n");

    // Iterate through events.
    uint64_t event_i = 0;
    while(sim_tree.Next()) {
        // Iterate through particles.
        for (int sh_it = 0; sh_it < idx.GetSize(); ++sh_it) {
            // Write to stdout.
            fprintf(
                f_parts, "%lu,%d,%d,%f,%lf,%lf,%lf,%f,%f,%f\n",
                event_i, idx[sh_it], pdg[sh_it], time[sh_it],
                vx[sh_it], vy[sh_it], vz[sh_it], px[sh_it], py[sh_it], pz[sh_it]
            );
        }

        ++event_i;
    }

    // Clean-up.
    fclose(f_parts);
    return 0;
}

/**
 * Extract hits from simulated and reconstructed ROOT files, associate, and wri-
 *     te them to a csv file.
 */
int write_hits(const char *f_sim, const char *f_rec) {
    // Initialize global dRICH dd4hep detector instance.
    init_dRICH();

    // Write hits and particles to output file.
    write_rec_hits( f_rec, f_sim, "csv/rec_hitmap.csv");
    write_sim_parts(f_sim, "csv/sim_hitmap.csv");

    return 0;
}

int input_handler() {
    write_hits(
        "/home/twig/code/eic/drich-dev/out/sim.edm4hep.root", // simu.
        "/home/twig/code/eic/drich-dev/out/rec.noise.edm4hep.root" // reco.
    );

    return 0;
}
