#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <charconv>
#include <thread>
#include <sstream>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>

using namespace std;
using boost::math::chi_squared;
using boost::math::normal;

// Pre-construct our standard normal so we don’t rebuild it each call
namespace {
    const normal NORMAL_DIST(0.0, 1.0);
}

// meta_analysis exactly as before
tuple<double,double,double,double,double,double,double,double,double>
meta_analysis(const vector<double>& beta, const vector<double>& se) {
    size_t n = beta.size();
    if (n == 0) return make_tuple(NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN);
    if (n == 1) {
        double mu = beta[0], se0 = se[0];
        double z  = se0>0 ? mu/se0 : 0.0;
        double p  = 2.0 * cdf(complement(NORMAL_DIST, fabs(z)));
        return make_tuple(p, se0, z, mu, 0.0, NAN, p, mu, z);
    }
    vector<double> var(n), w_fix(n);
    for (size_t i = 0; i < n; ++i) {
        var[i]   = se[i]*se[i];
        w_fix[i] = var[i]>0 ? 1.0/var[i] : 0.0;
    }
    double sum_w_fix = accumulate(w_fix.begin(), w_fix.end(), 0.0);
    if (sum_w_fix == 0.0) return make_tuple(NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN);

    double mu_fix = inner_product(w_fix.begin(), w_fix.end(), beta.begin(), 0.0)/sum_w_fix;
    double se_fix = sqrt(1.0/sum_w_fix);
    double z_fix  = se_fix>0 ? mu_fix/se_fix : 0.0;
    double p_fix  = 2.0 * cdf(complement(NORMAL_DIST, fabs(z_fix)));

    double Q = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double d = beta[i] - mu_fix;
        Q += w_fix[i] * d * d;
    }
    double df = double(n - 1);
    if (df <= 0 || isnan(Q))
        return make_tuple(NAN,NAN,NAN,NAN,NAN,NAN,p_fix,mu_fix,z_fix);

    double sum_w2 = 0.0;
    for (double w : w_fix) sum_w2 += w*w;
    double tau2 = max(0.0, (Q - df) / (sum_w_fix - sum_w2/sum_w_fix));

    vector<double> w_rand(n);
    for (size_t i = 0; i < n; ++i) {
        w_rand[i] = 1.0 / (var[i] + tau2);
    }
    double sum_w_rand = accumulate(w_rand.begin(), w_rand.end(), 0.0);
    if (sum_w_rand == 0.0)
        return make_tuple(NAN,NAN,NAN,NAN,NAN,NAN,p_fix,mu_fix,z_fix);

    double mu_rand = inner_product(w_rand.begin(), w_rand.end(), beta.begin(), 0.0)/sum_w_rand;
    double se_rand = sqrt(1.0/sum_w_rand);
    double z_rand  = se_rand>0 ? mu_rand/se_rand : 0.0;
    double p_rand  = 2.0 * cdf(complement(NORMAL_DIST, fabs(z_rand)));

    double I2     = (Q>df ? ((Q-df)/Q)*100.0 : 0.0);
    chi_squared chi2(df);
    double p_Q    = cdf(complement(chi2, Q));

    return make_tuple(
        p_rand, se_rand, z_rand, mu_rand,
        I2,     p_Q,     p_fix,  mu_fix, z_fix
    );
}

struct Accum {
    int chr, bp;
    vector<double> beta, se;
};

void process_files(const vector<string>& filenames, unsigned int n_threads) {
    if (filenames.empty()) {
        cerr << "Error: no input files provided.\n";
        return;
    }

    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    unordered_map<string, Accum> snp_map;
    static const array<string,5> req = { "SNP","CHR","BP","BETA","SE" };
    array<int,5> idx;

    // --- Sequential file reading & parsing ---
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
        vector<string> cols;
        for (size_t s = 0, p; (p = header.find('\t', s)) != string::npos; s = p+1)
            cols.emplace_back(header, s, p-s);
        cols.emplace_back(header, header.find_last_of('\t')+1);

        for (int i = 0; i < 5; ++i) {
            auto it = find(cols.begin(), cols.end(), req[i]);
            if (it == cols.end()) {
                cerr << "Error: missing column '" << req[i]
                     << "' in " << fn << "\n";
                return;
            }
            idx[i] = int(it - cols.begin());
        }

        string line;
        vector<string> fields;
        fields.reserve(cols.size());
        while (getline(in, line)) {
            if (line.empty()) continue;
            fields.clear();
            for (size_t s = 0, p; (p = line.find('\t', s)) != string::npos; s = p+1)
                fields.emplace_back(line, s, p-s);
            fields.emplace_back(line, line.find_last_of('\t')+1);
            if (fields.size() < cols.size()) continue;

            const string& snp = fields[idx[0]];
            int chr = strtol(fields[idx[1]].c_str(), nullptr, 10);
            int bp  = strtol(fields[idx[2]].c_str(), nullptr, 10);
            double b  = strtod(fields[idx[3]].c_str(), nullptr);
            double s0 = strtod(fields[idx[4]].c_str(), nullptr);

            auto& acc = snp_map[snp];
            if (acc.beta.empty()) {
                acc.chr = chr;
                acc.bp  = bp;
            }
            acc.beta.push_back(b);
            acc.se.push_back(s0);
        }
    }

    // --- Prepare parallel processing ---
    vector<string> keys;
    keys.reserve(snp_map.size());
    for (auto const& kv : snp_map) keys.push_back(kv.first);
    sort(keys.begin(), keys.end());

    cout << "SNP\tCHR\tBP\tN\tP\tSE\tBETA\tI2\tpQ\tBETA(FE)\tP(FE)\n";
    cout << scientific << setprecision(4);

    if (n_threads == 0) {
        n_threads = thread::hardware_concurrency();
        if (n_threads == 0) n_threads = 2;
    }
    size_t total = keys.size();
    size_t chunk = (total + n_threads - 1) / n_threads;

    vector<vector<string>> results(n_threads);
    vector<thread> workers;

    // --- Launch worker threads ---
    for (unsigned int t = 0; t < n_threads; ++t) {
        workers.emplace_back([&, t]() {
            ostringstream oss;
            oss << fixed;
            size_t start = t * chunk;
            size_t end   = min(total, start + chunk);
            for (size_t i = start; i < end; ++i) {
                auto const& snp = keys[i];
                auto const& A   = snp_map[snp];
                auto stats      = meta_analysis(A.beta, A.se);

                double p_rand  = get<0>(stats);
                double se_rand = get<1>(stats);
                double mu_rand = get<3>(stats);
                double I2      = get<4>(stats);
                double pQ      = get<5>(stats);
                double mu_fix  = get<7>(stats);
                double p_fix   = get<6>(stats);

                oss.str(""); oss.clear();
                oss << snp << '\t'
                    << A.chr << '\t'
                    << A.bp  << '\t'
                    << A.beta.size() << '\t'
                    << p_rand << '\t'
                    << se_rand << '\t'
                    << mu_rand << '\t'
                    << I2      << '\t'
                    << pQ      << '\t'
                    << mu_fix  << '\t'
                    << p_fix;
                results[t].push_back(oss.str());
            }
        });
    }

    // --- Join and print ---
    for (auto &th : workers) th.join();
    for (unsigned int t = 0; t < n_threads; ++t) {
        for (auto &line : results[t]) {
            cout << line << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0]
             << " <n_threads:0=auto> <input1> [<input2> ...]\n";
        return 1;
    }
    unsigned int n_threads = 0;
    try {
        n_threads = stoi(argv[1]);
    } catch (...) {
        cerr << "Error: first argument must be an integer (0 for auto).\n";
        return 1;
    }
    vector<string> files(argv+2, argv+argc);
    process_files(files, n_threads);
    return 0;
}
