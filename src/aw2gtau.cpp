#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstring>
#include <complex>
#include <random>
#include <functional>
#include <algorithm>

const double expcutoff = -708; // if N is smaller, exp(N) returns zero
const double one = 1;
const double two = 2;
const double pi = atan(1)*4;

int main() {
    std::string line;
    std::vector<std::string> vecString;
    std::vector<double> W;
    std::vector<double> Aw;
    char awfile[300];
    printf ("Input file name (spectral function):");
    if (scanf ("%299s",awfile) != 1) {printf("failed to read a value.\n");}

    std::fstream fin(awfile);
    if(fin.is_open()) {
        while (getline(fin, line)) {
                char *sArr = new char[line.length()+1];
                strcpy(sArr, line.c_str());
                char *pos = strtok(sArr," ");
                W.push_back(atof(pos));
                pos = strtok(NULL," ");
                Aw.push_back(atof(pos));
                delete[] sArr;
         }
    }
    else {
        std::cerr << "File could not be opened." << std::endl;
        return 1;
    }
    char outfile[300];
    printf ("Output file name:");
    if (scanf ("%299s",outfile) != 1) {
	printf("failed to read a value.\n");
    } else {
	std::ifstream output(outfile);
        if (output.good() == true) {
	    std::cerr << "WARNING: file already exists, will be rewritten.\n" << std::endl;
        }
    }	
    double beta;
    printf ("Inverse temperature (beta): ");
    if (scanf ("%lf",&beta) != 1) {printf("failed to read a value.\n");}
    int ntau;
    printf ("Number of imaginary time slices: ");
    if (scanf ("%i",&ntau) != 1) {printf("failed to read a value.\n");}
    double mu;
    printf ("Chemical potential: ");
    if (scanf ("%lf",&mu) != 1) {printf("failed to read a value.\n");}
    double noise_stdev;
    printf ("Standard deviation for noise in the real part of G (0 implies no noise): ");
    if (scanf ("%lf",&noise_stdev) != 1) {printf("failed to read a value.\n");}

    long nw = W.size();
    double step = ( W.back() - W.front() ) / nw;
    double tau = 0.;

    std::ofstream output;
    output.open(outfile, std::ofstream::out | std::ofstream::trunc);

//  random noise
    std::random_device seed;
    std::mt19937 engine(seed());
    std::normal_distribution<double> dist(1.0, noise_stdev);
    auto rnd = bind(dist, engine);
 
    std::vector<double> Noise(ntau);
    generate(Noise.begin(), Noise.end(), rnd);

    for (int z = 0; z != ntau; ++z) {
        std::vector<double> Kernel (nw);
        fill(Kernel.begin(), Kernel.end(), 0);
        for (int i = 0; i < nw; ++i) {
	    // avoiding NaN due to exp(N) for large negative N
            if ( W[i] * beta < expcutoff) Kernel[i] = 0.;
            else Kernel[i] = exp ( (mu - W[i]) * tau) / (one + exp ( (mu - W[i]) * beta) );
        }
        std::complex<double> s = 0;
        for (long i = 0; i < (nw - 1); ++i) {
            s += ( Aw[i] * Kernel[i] + Aw[i+1] * Kernel[i+1] ) / two;
        }
        output << tau << " " << step * s.real() * Noise[z] <<
                       " " << step * s.imag() << std::endl;
	tau += beta / (ntau - 1);
    }
    output.close();

    return 0;
}	
