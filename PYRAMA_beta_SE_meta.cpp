#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <charconv>         // for std::from_chars if you prefer
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>

using namespace std;
using boost::math::chi_squared;
using boost::math::normal;

// Pre‐construct our standard normal so we don’t rebuild it each call:
namespace {
    const normal NORMAL_DIST(0.0, 1.0);
}

// Computes both fixed‐ and random‐effects meta‐analysis; exactly as before,
// but renamed inputs and reusing the global NORMAL_DIST.
tuple<
    double, // p_random
    double, // se_random
    double, // z_random
    double, // mu_random
    double, // I2
    double, // p_Q
    double, // p_fixed
    double, // mu_fixed
    double  // z_fixed
> meta_analysis(const vector<double>& beta,
                const vector<double>& se)
{
    size_t n = beta.size();
    if (n == 0) {
        return make_tuple(NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN);
    }
    if (n == 1) {
        double mu = beta[0];
        double se0 = se[0];
        double z = se0 > 0 ? mu / se0 : 0.0;
        double p = 2.0 * cdf(complement(NORMAL_DIST, fabs(z)));
        return make_tuple(p, se0, z, mu, 0.0, NAN, p, mu, z);
    }

    // Precompute variances and fixed weights
    vector<double> var(n), w_fix(n);
    for (size_t i = 0; i < n; ++i) {
        var[i]    = se[i]*se[i];
        w_fix[i]  = var[i] > 0.0 ? 1.0/var[i] : 0.0;
    }
    double sum_w_fix = accumulate(w_fix.begin(), w_fix.end(), 0.0);
    if (sum_w_fix == 0.0) {
        return make_tuple(NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN);
    }

    // Fixed‐effect
    double mu_fix = inner_product(w_fix.begin(), w_fix.end(), beta.begin(), 0.0) / sum_w_fix;
    double se_fix = sqrt(1.0 / sum_w_fix);
    double z_fix  = se_fix > 0 ? mu_fix / se_fix : 0.0;
    double p_fix  = 2.0 * cdf(complement(NORMAL_DIST, fabs(z_fix)));

    // Cochran's Q
    double Q = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double d = beta[i] - mu_fix;
        Q += w_fix[i] * d * d;
    }
    double df = double(n - 1);
    if (df <= 0 || isnan(Q)) {
        return make_tuple(NAN,NAN,NAN,NAN,NAN,NAN,p_fix,mu_fix,z_fix);
    }

    // DerSimonian–Laird tau^2
    double sum_w2 = 0.0;
    for (double w : w_fix) sum_w2 += w*w;
    double tau2 = max(0.0, (Q - df) / (sum_w_fix - sum_w2/sum_w_fix));

    // Random‐effects weights
    vector<double> w_rand(n);
    for (size_t i = 0; i < n; ++i) {
        w_rand[i] = 1.0 / (var[i] + tau2);
    }
    double sum_w_rand = accumulate(w_rand.begin(), w_rand.end(), 0.0);
    if (sum_w_rand == 0.0) {
        return make_tuple(NAN,NAN,NAN,NAN,NAN,NAN,p_fix,mu_fix,z_fix);
    }

    double mu_rand = inner_product(w_rand.begin(), w_rand.end(), beta.begin(), 0.0) / sum_w_rand;
    double se_rand = sqrt(1.0 / sum_w_rand);
    double z_rand  = se_rand > 0 ? mu_rand / se_rand : 0.0;
    double p_rand  = 2.0 * cdf(complement(NORMAL_DIST, fabs(z_rand)));

    double I2    = (Q > df ? ((Q - df)/Q) * 100.0 : 0.0);
    chi_squared chi2(df);
    double p_Q   = cdf(complement(chi2, Q));

    return make_tuple(
        p_rand, se_rand, z_rand, mu_rand,
        I2, p_Q, p_fix, mu_fix, z_fix
    );
}

// Holds the accumulating data for each SNP
struct Accum {
    int    chr, bp;
    vector<double> beta, se;
};

void process_files(const vector<string>& filenames) {
    if (filenames.empty()) {
        cerr << "Error: no input files provided.\n";
        return;
    }

    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Our per‐SNP accumulator:
    unordered_map<string, Accum> snp_map;
    // Reserve to avoid rehashing if you have a rough idea of SNP count:
    // snp_map.reserve(1 << 20);

    // Required column names:
    static const array<string,5> req = { "SNP","CHR","BP","BETA","SE" };
    array<int,5> idx;  // positions of those columns in the current file

    for (auto const& fn : filenames) {
        ifstream in(fn);
        if (!in.is_open()) {
            cerr << "Error: cannot open '" << fn << "'\n";
            return;
        }

        string header;
        getline(in, header);
        if (header.empty()) {
            cerr << "Error: '" << fn << "' is empty\n";
            return;
        }

        // Split header on tabs (manual, faster than istringstream)
        vector<string> cols;
        cols.reserve(16);
        size_t start = 0, pos;
        while ((pos = header.find('\t', start)) != string::npos) {
            cols.emplace_back(header, start, pos - start);
            start = pos + 1;
        }
        cols.emplace_back(header, start);

        // Find our required columns
        for (int i = 0; i < 5; ++i) {
            auto it = find(cols.begin(), cols.end(), req[i]);
            if (it == cols.end()) {
                cerr << "Error: missing column '" << req[i]
                     << "' in " << fn << "\n";
                return;
            }
            idx[i] = int(it - cols.begin());
        }

        // Parse each line
        string line;
        vector<string> fields;
        fields.reserve(cols.size());
        while (getline(in, line)) {
            if (line.empty()) continue;
            // split on tabs:
            fields.clear();
            size_t s = 0, p;
            while ((p = line.find('\t', s)) != string::npos) {
                fields.emplace_back(line, s, p - s);
                s = p + 1;
            }
            fields.emplace_back(line, s);

            if ((int)fields.size() < (int)cols.size()) continue;

            const string& snp = fields[idx[0]];
            // fast parse chr, bp, beta, se:
            int    chr = strtol(fields[idx[1]].c_str(), nullptr, 10);
            int    bp  = strtol(fields[idx[2]].c_str(), nullptr, 10);
            double b   = strtod(fields[idx[3]].c_str(), nullptr);
            double se  = strtod(fields[idx[4]].c_str(), nullptr);

            auto& acc = snp_map[snp];
            if (acc.beta.empty()) {
                acc.chr = chr;
                acc.bp  = bp;
            }
            acc.beta.push_back(b);
            acc.se.push_back(se);
        }
    }

    // Output header:
    cout << "SNP\tCHR\tBP\tN\t"
         << "P\tSE\tBETA\tI2\tpQ\t"
         << "BETA(FE)\tP(FE)\n";

    // To print sorted by SNP name:
    vector<string> keys;
    keys.reserve(snp_map.size());
    for (auto const& kv : snp_map) keys.push_back(kv.first);
    sort(keys.begin(), keys.end());

    // Compute & print
    cout << scientific << setprecision(4);
    for (auto const& snp : keys) {
        auto const& A = snp_map[snp];
        auto stats = meta_analysis(A.beta, A.se);

        double p_rand  = get<0>(stats);
        double se_rand = get<1>(stats);
        double mu_rand = get<3>(stats);
        double I2      = get<4>(stats);
        double pQ      = get<5>(stats);
        double p_fix   = get<6>(stats);
        double mu_fix  = get<7>(stats);

        cout << snp << '\t'
             << A.chr << '\t'
             << A.bp  << '\t'
             << int(A.beta.size()) << '\t'
             << p_rand << '\t'
             << defaultfloat << se_rand << '\t'
             << mu_rand << '\t'
             << I2      << '\t'
             << scientific << pQ << '\t'
             << defaultfloat << mu_fix << '\t'
             << scientific << p_fix << "\n";
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0]
             << " <input1> [<input2> ...]\n";
        return 1;
    }
    vector<string> files(argv+1, argv+argc);
    process_files(files);
    return 0;
}
