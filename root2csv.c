// C.
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// C++.
#include <list>
#include <map>

// ROOT.
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"

// dd4hep.
#include <DD4hep/Detector.h>
#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>

/* TODO LIST:
 *   * Implement a feature to process many files.
 *    ---> To do this, we'll have to think about how to associate each simula-
 *         ted with each reconstructed file.
 *   * Upload the ipython script to generate plots in the global coordinate sys-
 *         tem
 */

/** dRICH system number for the cell ID. */
#define DRICH_SYSTEM 120

/** String handling constants. */
#define CELLMAP_LOC "csv/cellmap.csv"
#define BUFFERSIZE 256

/** Variables to go from global to local coordinates. */
#define DILATION 4.3
#define OFFSET   0.5

const char *USAGE_MESSAGE =
    "\nUsage: root2csv [-hs:r:S:R:]"
    "\n * -h           : show this message and exit."
    "\n * -s f_in_sim  : input simulation file."
    "\n * -r f_in_rec  : input reconstructed file."
    "\n * -S f_out_sim : output sim csv file. Default is `csv/sim_parts.csv`."
    "\n * -R f_out_rec : output rec csv file. Default is `csv/rec_hits.csv`."
    "\n\n    Produce two input csv files from sim and rec root files.\n\n";

char DEF_SIMNAME[] = "csv/sim_parts.csv";
char DEF_RECNAME[] = "csv/rec_hits.csv";

/** Struct to store one reconstructed hit. */
struct rec_hit {
    uint64_t event;
    uint64_t sector;
    uint64_t x;
    uint64_t y;
    int32_t time;
    int32_t charge;
    int32_t pindex;
};

/** Global instances of the dd4hep dRICH detector and its bitfield coder. */
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
static double global_to_local(double d) {return DILATION*d + OFFSET;}

/** Check if d is in list l. */
static int32_t is_in_list(std::list<int32_t> l, int32_t d) {
    return std::find(l.begin(), l.end(), d) == l.end() ? 0 : 1;
}

/**
 * Initialize instances of dRICH detector and bit field coder from dd4hep. This
 *     function assumes that we are running over the eic-shell with the dRICH
 *     environ.sh script loaded.
 */
static int32_t init_dRICH() {
    dd4hep::Detector *det = &(dd4hep::Detector::getInstance());
    det->fromXML("/opt/detector/epic-23.10.0/share/epic/epic.xml");
    g_dRICH = det->detector("DRICH");
    g_readout_coder = det->readout("DRICHHits").idSpec().decoder();

    return 0;
}

/** Get minimum pixel x,y positions and store in x_min and y_min. */
static int32_t get_min_pix(int32_t *x_min, int32_t *y_min) {
    // Initialize to "infinity".
    *x_min = INT_MAX;
    *y_min = INT_MAX;

    // Iterate through dRICH sector 0 sensors.
    for (auto const &[d_name, d_sensor] : g_dRICH.children()) {
        if (d_name.find("sensor_de_sec0") == std::string::npos) continue;
        const auto *det_pars =
            d_sensor.extension<dd4hep::rec::VariantParameters>(true);

        // Transform pixel position to local coordinates.
        double x = global_to_local(det_pars->get<double>("pos_x"));
        double y = global_to_local(det_pars->get<double>("pos_y"));

        // If lower, set new minimum.
        if (x < (double) *x_min) *x_min = (int32_t) x;
        if (y < (double) *y_min) *y_min = (int32_t) y;
    }

    return 0;
}

/**
 * Print a cell map to a csv file under the name CELLMAP_LOC. The cell map index
 *     is the cell ID of each photomultiplier pixel in the sector 0 of dRICH.
 *     The stored data is the position of the pixel in the local coordinate sys-
 *     tem of each dRICH sector.
 */
static int32_t create_cellmap() {
    // Create cellmap file.
    FILE *f_map = fopen(CELLMAP_LOC, "w");

    // To define the 0,0 point of the matrix, find the minimum pixel positions.
    int32_t pixel_x_min, pixel_y_min;
    get_min_pix(&pixel_x_min, &pixel_y_min);

    // Write header.
    fprintf(f_map, "cell_id,x,y\n");

    // Iterate through the dRICH sub-detectors.
    for (auto const &[d_name, d_sensor] : g_dRICH.children()) {
        // Since the coordinate system is the same for all sectors, we only care
        //     about sector 0 here.
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
        uint64_t idx_x = (uint64_t) (((int32_t) pixel_x) - pixel_x_min);
        uint64_t idx_y = (uint64_t) (((int32_t) pixel_y) - pixel_y_min);

        // Walk through the set of 64 pixels.
        for (int32_t xi = 0; xi < 8; ++xi) {
            for (int32_t yi = 0; yi < 8; ++yi) {
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
static int32_t read_cellmap(
    std::map<uint64_t, std::pair<uint64_t, uint64_t>> *cellmap
) {
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
        for (int32_t val_i = 0; val != NULL; ++val_i) {
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
static int32_t write_rec_hits(
    const char *in_rec, const char *in_sim, const char *out
) {
    // Read cell map from CELLMAP_LOC into an std::map.
    std::map<uint64_t, std::pair<uint64_t, uint64_t>> cellmap;
    read_cellmap(&cellmap);

    // Open trees.
    TTreeReader rec_tree((TTree *) (new TFile(in_rec))->Get("events"));
    TTreeReader sim_tree((TTree *) (new TFile(in_sim))->Get("events"));

    // Associate TTreeReaderArrays with relevant data from trees.
    TTreeReaderArray<uint64_t> rec_cell_id(rec_tree, "DRICHRawHits.cellID");
    TTreeReaderArray<int32_t>  rec_charge (rec_tree, "DRICHRawHits.charge");
    TTreeReaderArray<int32_t>  rec_time   (rec_tree, "DRICHRawHits.timeStamp");

    TTreeReaderArray<uint64_t> sim_cell_id(sim_tree, "DRICHHits.cellID");
    TTreeReaderArray<int32_t>  sim_index  (sim_tree, "DRICHHits#0.index");

    // Set TTreeReaders to first entry.
    sim_tree.SetEntry(-1);
    rec_tree.SetEntry(-1);

    // Open hitmap output file and print header.
    FILE *f_map = fopen(out, "w");
    fprintf(f_map, "event,time,sector,x,y,charge,pindex\n");

    // Iterate through events.
    uint64_t event_i = 0;
    while(sim_tree.Next() && rec_tree.Next()) {
        // List to store all reconstructed hits in event.
        std::list<rec_hit> rec_hit_list;

        // Iterate through rec and sim hits.
        for (int32_t rh_it = 0; rh_it < rec_cell_id.GetSize(); ++rh_it) {
            for (int32_t sh_it = 0; sh_it < sim_cell_id.GetSize(); ++sh_it) {
                if (sim_cell_id[sh_it] != rec_cell_id[rh_it]) continue;
                uint64_t cell_id = rec_cell_id[rh_it];

                // Store sector and check hit position from cellmap.
                uint64_t sector = g_readout_coder->get(cell_id, "sector");
                g_readout_coder->set(cell_id, "sector", 0);

                if (!cellmap.count(cell_id)) {
                    fprintf(stderr, "%lu not in cellmap. Exiting.\n", cell_id);
                    return 1;
                }

                // Add hit to event list.
                rec_hit_list.push_back(rec_hit{
                    .event  = event_i,
                    .sector = sector,
                    .x      = cellmap[cell_id].first,
                    .y      = cellmap[cell_id].second,
                    .time   = rec_time[rh_it],
                    .charge = rec_charge[rh_it],
                    .pindex = sim_index[sh_it]
                });
            }
        }

        // Sort list based on timestamps.
        rec_hit_list.sort([](const rec_hit &a, const rec_hit &b) {
            return a.time < b.time;
        });

        // Print list to f_map.
        for (auto it = rec_hit_list.begin(); it != rec_hit_list.end(); ++it) {
            fprintf(
                f_map, "%lu,%d,%lu,%lu,%lu,%d,%d\n", it->event, it->time,
                it->sector, it->x, it->y, it->charge, it->pindex
            );
        }

        ++event_i;
    }

    // Clean-up.
    fclose(f_map);
    return 0;
}

/** Record all particle indices in in_csv into std::map. */
static int32_t record_indices(
    const char *in_csv, std::map<uint64_t, std::list<int32_t>> *pindices
) {
    FILE *f_hits = fopen(in_csv, "r");

    // Iterate through hitmap csv.
    uint64_t event = UINT_MAX;

    char buffer[BUFFERSIZE];
    fgets(buffer, BUFFERSIZE, f_hits); // Ignore header.
    while (fgets(buffer, BUFFERSIZE, f_hits)) {
        // First val is event, 5 are unrelated, and final val is pindex.
        char *val = strtok(buffer, ",");
        if (event != strtoull(val, NULL, 0)) {
            event = strtoull(val, NULL, 0);
            (*pindices)[event] = {};
        }

        // Move to last value.
        for (int32_t val_i = 0; val_i < 6; ++val_i) val = strtok(NULL, ",");

        // Write pindex to list.
        (*pindices)[event].push_back(strtol(val, NULL, 0));
    }

    // Clean-up
    fclose(f_hits);
    return 0;
}

/** Write simulated particles into a csv file. */
static int32_t write_sim_parts(
    const char *in_sim, const char *out,
    std::map<uint64_t, std::list<int32_t>> pindices
) {
    // Open sim tree.
    TTreeReader sim_tree((TTree *) (new TFile(in_sim))->Get("events"));

    // Associate TTreeReaderArrays with relevant data from tree.
    TTreeReaderArray<int32_t> idx (sim_tree, "DRICHHits#0.index");
    TTreeReaderArray<int32_t> pdg (sim_tree, "MCParticles.PDG");
    TTreeReaderArray<float>   time(sim_tree, "MCParticles.time");
    TTreeReaderArray<double>  vx  (sim_tree, "MCParticles.vertex.x");
    TTreeReaderArray<double>  vy  (sim_tree, "MCParticles.vertex.y");
    TTreeReaderArray<double>  vz  (sim_tree, "MCParticles.vertex.z");
    TTreeReaderArray<float>   px  (sim_tree, "MCParticles.momentum.x");
    TTreeReaderArray<float>   py  (sim_tree, "MCParticles.momentum.y");
    TTreeReaderArray<float>   pz  (sim_tree, "MCParticles.momentum.z");

    // Rewind TTreeReader.
    sim_tree.SetEntry(-1);

    // Open output file and print header.
    FILE *f_parts = fopen(out, "w");
    fprintf(f_parts, "event,pindex,pdg,time,vx,vy,vz,px,py,pz\n");

    // Iterate through events.
    uint64_t event_i = 0;
    while(sim_tree.Next()) {
        // Iterate through particles.
        for (int32_t sh_it = 0; sh_it < idx.GetSize(); ++sh_it) {
            // Check if index is in the std::map.
            if (!is_in_list(pindices[event_i], idx[sh_it])) continue;

            // Write to stdout.
            fprintf(
                f_parts, "%lu,%d,%d,%f,%lf,%lf,%lf,%f,%f,%f\n", event_i,
                idx[sh_it], pdg[sh_it], time[sh_it], vx[sh_it], vy[sh_it],
                vz[sh_it], px[sh_it], py[sh_it], pz[sh_it]
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
static int32_t write_hits(
    char *f_in_sim, char *f_in_rec, char *f_out_sim, char *f_out_rec
) {
    // Initialize global dRICH dd4hep detector instance.
    init_dRICH();

    // Write reconstructed hits to csv.
    write_rec_hits(f_in_rec, f_in_sim, f_out_rec);

    // Record all pindices associated to reconstructed hits.
    std::map<uint64_t, std::list<int32_t>> pindex_lists;
    record_indices(f_out_rec, &pindex_lists);

    // Write simulated particles whose hits were reconstructed.
    write_sim_parts(f_in_sim, f_out_sim, pindex_lists);

    return 0;
}

static int32_t copy_str(char *src, char **tgt) {
    *tgt = (char *) (malloc(strlen(src) + 1));
    strcpy(*tgt, src);
    return 0;
}

/** Handle arguments using optarg. */
static int32_t handle_args(
    int32_t argc, char **argv,
    char **f_in_sim, char **f_in_rec, char **f_out_sim, char **f_out_rec
) {
    int32_t opt;
    while ((opt = getopt(argc, argv, "-hs:r:S:R:")) != -1) {
        switch (opt) {
            case 'h':
                return 2;
            case 's':
                copy_str(optarg, f_in_sim);
                break;
            case 'r':
                copy_str(optarg, f_in_rec);
                break;
            case 'S':
                copy_str(optarg, f_out_sim);
                break;
            case 'R':
                copy_str(optarg, f_out_rec);
                break;
            default:
                return 1;
        }
    }

    // Check that input files were given.
    if (*f_in_sim == NULL || *f_in_rec == NULL) {
        fprintf(stderr, "Error: Please provide input files to convert.\n");
        return 1;
    }

    // Fill output files if empty.
    if (*f_out_sim == NULL) copy_str(DEF_SIMNAME, f_out_sim);
    if (*f_out_rec == NULL) copy_str(DEF_RECNAME, f_out_rec);

    return 0;
}

/** Print usage and return error. */
static int32_t usage() {
    printf("%s", USAGE_MESSAGE);
    return 1;
}

/** Entry point of the program. */
int32_t main(int32_t argc, char **argv) {
    // Handle input.
    char *f_in_sim  = NULL;
    char *f_in_rec  = NULL;
    char *f_out_sim = NULL;
    char *f_out_rec = NULL;

    if (handle_args(argc, argv, &f_in_sim, &f_in_rec, &f_out_sim, &f_out_rec)) {
        return usage();
    }

    // Run.
    write_hits(f_in_sim, f_in_rec, f_out_sim, f_out_rec);

    // Free up memory.
    if (f_in_sim  != NULL) free(f_in_sim);
    if (f_in_rec  != NULL) free(f_in_rec);
    if (f_out_sim != NULL) free(f_out_sim);
    if (f_out_rec != NULL) free(f_out_rec);

    return 0;
}
