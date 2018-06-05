#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstring>
#include <complex>
#include <stdlib.h>

const double two = 2;
const double pi = atan(1)*4;
const std::complex<double> j(0, 1);

// integrand
std::complex<double> f(double  iw,
                          int  w,
                       double  mu,
           std::vector<double> *W,
           std::vector<double> *Aw) {
    return (*Aw)[w] / (j * iw - (*W)[w] + mu);}

int main(int argc, char* argv[]) {
    if (argc > 1) {
        std::cerr << "An interactive program which transforms a given A(w) " 
             << " [a two-column file] into G(iw)" << std::endl;
        return 1;
    }
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
    int niw;
    printf ("Number of Matsubara frequencies: ");
    if (scanf ("%i",&niw) != 1) {printf("failed to read a value.\n");}
    double mu;
    printf ("Chemical potential: ");
    if (scanf ("%lf",&mu) != 1) {printf("failed to read a value.\n");}

    std::ofstream output;
    output.open (outfile, std::ofstream::out | std::ofstream::trunc);
    int nw = W.size();
    double step = ( W.back() - W.front() ) / nw;

    output << "#Column 1: iw" << std::endl 
           << "#Column 2: Re[G(iw)]" << std::endl
           << "#Column 3: Im[G(iw)]" << std::endl;
    for (int x = 0; x != niw; ++x) {
        double iw = (2*x + 1 - niw) * pi / beta;
    	std::complex<double> s = 0;
        for (int i = 0; i < (nw - 1); ++i)
            s += ( f(iw, i+1, mu, &W, &Aw) + f(iw, i, mu, &W, &Aw) ) / two;
        output << iw << " " << step * s.real() <<
                        " " << step * s.imag() << std::endl;
    }
    output.close();

    return 0;
}
