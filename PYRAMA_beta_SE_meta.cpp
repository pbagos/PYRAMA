#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <string>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip> // For output precision
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>

using namespace std;
using boost::math::chi_squared;
using boost::math::normal;

// Define a struct to hold SNP data
struct SNPData {
    string snp;
    int chr;
    int bp;
    double beta;
    double se;
};

// Function to calculate meta-analysis statistics.
// Returns a 9‐element tuple (we keep computing Z internally, but will not print it):
// (p_random, se_random, z_random, mu_random, I2, p_Q, p_fixed, mu_fixed, z_fixed)
tuple<double, double, double, double, double, double, double, double, double>
altmeta(const vector<double>& y1, const vector<double>& s2) {
    size_t n = y1.size();

    // If no studies, return all NaN
    if (n == 0) {
        return make_tuple(
            NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN
        );
    }

    // If exactly one study, random and fixed coincide
    if (n == 1) {
        double mu = y1[0];
        double se = s2[0];
        double z = (se > 0 ? mu / se : 0.0);
        normal normal_dist(0.0, 1.0);
        double p = 2.0 * cdf(complement(normal_dist, fabs(z)));

        return make_tuple(
            p,     // p_random
            se,    // se_random
            z,     // z_random
            mu,    // mu_random  → BETA
            0.0,   // I_squared
            NAN,   // p_Q (undefined for single study)
            p,     // p_fixed
            mu,    // mu_fixed  → BETA(FE)
            z      // z_fixed
        );
    }

    // More than one study: compute variances and fixed weights
    vector<double> variances(n), weights_fixed(n), weights_random(n);
    for (size_t i = 0; i < n; ++i) {
        variances[i] = s2[i] * s2[i];
        weights_fixed[i] = (variances[i] > 0 ? 1.0 / variances[i] : 0.0);
    }

    double weights_sum_fixed = accumulate(
        weights_fixed.begin(), weights_fixed.end(), 0.0
    );
    if (weights_sum_fixed == 0.0) {
        // All variances zero or invalid
        return make_tuple(
            NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN
        );
    }

    // Fixed‐effect combined mean (mu_fixed) and its SE
    double mu_fixed = inner_product(
        weights_fixed.begin(), weights_fixed.end(),
        y1.begin(), 0.0
    ) / weights_sum_fixed;
    double se_fixed = sqrt(1.0 / weights_sum_fixed);
    double z_fixed = (se_fixed > 0 ? mu_fixed / se_fixed : 0.0);
    normal normal_dist(0.0, 1.0);
    double p_fixed = 2.0 * cdf(complement(normal_dist, fabs(z_fixed)));

    // Cochran's Q for heterogeneity
    double Q = 0.0;
    for (size_t i = 0; i < n; ++i) {
        Q += weights_fixed[i] * pow(y1[i] - mu_fixed, 2);
    }

    double df = static_cast<double>(n - 1);
    if (df <= 0 || isnan(Q)) {
        // Not enough degrees of freedom
        return make_tuple(
            NAN,       // p_random
            NAN,       // se_random
            NAN,       // z_random
            NAN,       // mu_random
            NAN,       // I_squared
            NAN,       // p_Q
            p_fixed,   // p_fixed
            mu_fixed,  // mu_fixed → BETA(FE)
            z_fixed    // z_fixed
        );
    }

    // Between‐study variance tau^2 (DerSimonian‐Laird)
    double sum_w2 = 0.0;
    for (double w : weights_fixed) {
        sum_w2 += w * w;
    }
    double tau_squared = max(
        0.0,
        (Q - df) / (weights_sum_fixed - (sum_w2 / weights_sum_fixed))
    );

    // Compute random‐effects weights
    for (size_t i = 0; i < n; ++i) {
        weights_random[i] = 1.0 / (variances[i] + tau_squared);
    }

    double weights_sum_random = accumulate(
        weights_random.begin(), weights_random.end(), 0.0
    );
    if (weights_sum_random == 0.0) {
        // Cannot compute random‐effects
        return make_tuple(
            NAN,       // p_random
            NAN,       // se_random
            NAN,       // z_random
            NAN,       // mu_random
            NAN,       // I_squared
            NAN,       // p_Q
            p_fixed,   // p_fixed
            mu_fixed,  // mu_fixed → BETA(FE)
            z_fixed    // z_fixed
        );
    }

    // Random‐effects combined mean (mu_random) and its SE
    double mu_random = inner_product(
        weights_random.begin(), weights_random.end(),
        y1.begin(), 0.0
    ) / weights_sum_random;
    double se_random = sqrt(1.0 / weights_sum_random);
    double z_random = (se_random > 0 ? mu_random / se_random : 0.0);
    double p_random = 2.0 * cdf(complement(normal_dist, fabs(z_random)));

    // Heterogeneity I^2 and p_Q
    double I_squared = (Q > df) ? ((Q - df) / Q) * 100.0 : 0.0;
    chi_squared chi_dist(df);
    double p_Q = cdf(complement(chi_dist, Q));

    return make_tuple(
        p_random,  // p_random
        se_random, // se_random
        z_random,  // z_random
        mu_random, // mu_random → BETA
        I_squared, // I^2
        p_Q,       // heterogeneity p‐value
        p_fixed,   // p_fixed
        mu_fixed,  // mu_fixed → BETA(FE)
        z_fixed    // z_fixed
    );
}

// Main function to read data and process meta-analysis
void processFile(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return;
    }

    string line;
    getline(file, line); // Read the header

    // Parse header to find indices of required columns
    istringstream headerStream(line);
    vector<string> headerColumns;
    string column;
    while (headerStream >> column) {
        headerColumns.push_back(column);
    }

    // Map column names to indices
    unordered_map<string, int> columnIndices;
    for (size_t i = 0; i < headerColumns.size(); ++i) {
        columnIndices[headerColumns[i]] = i;
    }

    // Ensure required columns exist
    vector<string> requiredColumns = { "SNP", "CHR", "BP", "BETA", "SE" };
    for (const auto& col : requiredColumns) {
        if (columnIndices.find(col) == columnIndices.end()) {
            cerr << "Error: Required column " << col << " not found in the input file." << endl;
            return;
        }
    }

    vector<SNPData> snpData;

    // Read and store each data line (tab-delimited)
    while (getline(file, line)) {
        istringstream lineStream(line);
        vector<string> rowData;
        string value;

        while (getline(lineStream, value, '\t')) {
            rowData.push_back(value);
        }

        // Basic sanity check: skip incomplete lines
        if (rowData.size() < headerColumns.size()) continue;

        SNPData snp;
        snp.snp  = rowData[columnIndices["SNP"]];
        snp.chr  = stoi(rowData[columnIndices["CHR"]]);
        snp.bp   = stoi(rowData[columnIndices["BP"]]);
        snp.beta = stod(rowData[columnIndices["BETA"]]);
        snp.se   = stod(rowData[columnIndices["SE"]]);

        snpData.push_back(snp);
    }

    // Group betas and ses by SNP
    unordered_map<string, vector<double>> beta_map, se_map;
    unordered_map<string, int> chr_map, bp_map;
    for (const auto& snp : snpData) {
        beta_map[snp.snp].push_back(snp.beta);
        se_map[snp.snp].push_back(snp.se);
        chr_map[snp.snp] = snp.chr;
        bp_map[snp.snp]  = snp.bp;
    }

    // ------------- CHANGED HEADER -------------
    // Print new header WITH an N column:
    // SNP  CHR  BP  N  P  SE  BETA  I2  pQ  BETA(FE)  P(FE)
    cout << "SNP\tCHR\tBP\tN\t"
            "P\tSE\tBETA\tI2\tpQ\t"
            "BETA(FE)\tP(FE)"
         << endl;

    // Loop over each SNP
    for (const auto& kv : beta_map) {
        const string& snp_name = kv.first;
        const auto& beta_list = kv.second;
        const auto& se_list   = se_map.at(snp_name);
        int chr = chr_map[snp_name];
        int bp  = bp_map[snp_name];

        // ------------- COMPUTE N -------------
        int N = static_cast<int>(beta_list.size());

        auto result = altmeta(beta_list, se_list);
        double p_random   = get<0>(result);
        double se_random  = get<1>(result);
        double mu_random  = get<3>(result); // → BETA (random)
        double I2         = get<4>(result);
        double p_Q        = get<5>(result);
        double p_fixed    = get<6>(result);
        double mu_fixed   = get<7>(result); // → BETA(FE)

        // Output fields in the order:
        // snp, chr, bp, N,
        // p_random, se_random, mu_random, I2, p_Q,
        // mu_fixed (BETA(FE)), p_fixed (P(FE))
        cout << snp_name << "\t"
             << chr       << "\t"
             << bp        << "\t"
             << N         << "\t";  // <-- print N here

        // Random‐effect p (scientific, precision 4)
        cout << scientific << setprecision(4) << p_random << "\t";

        // Then default formatting for the rest:
        cout << defaultfloat
             << se_random  << "\t"
             << mu_random  << "\t"
             << I2         << "\t";

        // Heterogeneity p_Q (scientific, precision 4)
        cout << scientific << setprecision(4) << p_Q << "\t";

        // Fixed‐effect columns (contiguous):
        // BETA(FE), P(FE)
        cout << defaultfloat
             << mu_fixed << "\t";

        cout << scientific << setprecision(4) << p_fixed << endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        return 1;
    }

    string filename = argv[1];
    processFile(filename);
    return 0;
}
