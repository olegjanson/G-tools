#include "../include/hk2gwfw.h"

const double pi = atan(1)*4;
const std::complex<double> j(0, 1);

void downfold(Eigen::MatrixXcd *M,
              Eigen::MatrixXcd *Md,
                           int  blocksize){
//  downfolding of a matrix: diagonal blocks of the size `blocksize`
//  are copied from M to Md. If `blocksize`=1, only the diagonal is copied.
    int nblocks = (*M).rows() / blocksize;
    for (int bi = 0; bi < nblocks; ++bi){
        (*Md).block(bi*blocksize,
                    bi*blocksize,
                    blocksize,
                    blocksize) = ((*M).block(bi*blocksize,
                                             bi*blocksize,
                                             blocksize,
                                             blocksize));
    }
    return;
}

void blockwise_inverse(Eigen::MatrixXcd *M,
                       Eigen::MatrixXcd *Minv,
                                    int  blocksize){
//  a blockwise inverse of the matrix `m`: returns a block-diagonal matrix `m'`
//  whose blocks of the size `blocksize` are inverses of the respective blocks
//  of `m`, all other matrix elements are zero. If `blocksize=1`, `m'` is a
//  diagonal matrix with the matrix elements given by (1 / m_ii).
    int nblocks = (*M).rows() / blocksize;
    for (int bi = 0; bi < nblocks; ++bi){
        (*Minv).block(bi*blocksize,
                      bi*blocksize,
                      blocksize,
                      blocksize) = (((*M).block(bi*blocksize,
                                                bi*blocksize,
                                                blocksize,
                                                blocksize)
                                    ).inverse());
    }
    return;
}

int main (int argc, char* argv[]) {

    std::cout << "This interactive code evaluates the (non-interacting) ";
    std::cout << "lattice Green's function G(w), the spectral function A(w), ";
    std::cout << "as well as hybridization function F(w) for a multi-orbital ";
    std::cout << "Hamiltonian given in the *.hk format. The frequencies can ";
    std::cout << "be real or Matsubara." << std::endl;
    std::cout << std::endl << "Please report bugs to Oleg Janson ";
    std::cout << "<oleg" << "janson" << "@gmail.com>." << std::endl << std::endl;

//  INPUT
    char infile [300];
    printf ("Input file name (*.hk file):");
    if (scanf ("%299s",infile) != 1) {
        printf("failed to read a value.\n");
    } else {
        std::ifstream hkfile(infile);
        if (hkfile.good() != true) {
	    std::cerr << "Failed to find the file." << std::endl;
	    return 1;
        }
    }

    char outfile [300];
    printf ("Output file name:");
    if (scanf ("%299s",outfile) != 1) {
        printf("failed to read a value.\n");
    } else {
        std::ifstream output(outfile);
        if (output.good() == true) {
	    std::cerr << "WARNING: file already exists, will be rewritten.\n" << std::endl;
        }
    }	

    bool iw = false;
    double wmin, wmax, beta, mureal;
    int nw;
    printf("Lowest frequency (set 0 for the Matsubara axis): ");
    if (scanf("%lf",&wmin) != 1) {printf("failed to read a value.\n");}
    printf ("Highest frequency (set 0 for the Matsubara axis): ");
    if (scanf("%lf",&wmax) != 1) {printf("failed to read a value.\n");}
    if ((wmin == 0) && (wmax == 0)) {
        printf ("You have chosen the Matsubara axis.  Specify the reverse temperature (beta): ");
        if (scanf("%lf",&beta) != 1) {printf("failed to read a value.\n");}
        iw = true;
    }
    printf ("Number of frequencies (real or positive Matsubara): ");
    if (scanf("%i",&nw) != 1) {printf("failed to read a value.\n");}
    printf ("The chemical potential: ");
    if (scanf("%lf",&mureal) != 1) {printf("failed to read a value.\n");}
    std::complex<double> mu = mureal;
    std::complex<double> idelta = 3. * (wmax - wmin) * j * (1. / nw); 

    char answer;
    printf("Which matrix elements to print: diagonal (d), nondiagonal (n), or both (b)?");
    if (scanf (" %c",&answer) != 1) {printf("failed to read a value.\n");}
    bool printDiag = false;
    bool printNonDiag = false;
    if (answer == 'd') {
 	printDiag = true;
    } else if (answer == 'n') {
 	printNonDiag = true;
    } else if (answer == 'b') {
        printDiag = true;
        printNonDiag = true;
    } else {
        std::cerr << "Did not get it!" << std::endl;
        return 1;
    }

    printf("What to print: hybridization function (d), Green's function (g) or spectral function (a)?");
    if (scanf (" %c",&answer) != 1) {printf("failed to read a value.\n");}
    bool printDelta = false;
    bool doDownfolding = false;
    int blocksize;
    bool printG = false;
    bool printA = false;
    if (answer == 'd') {
        printDelta = true;
    } else if (answer == 'g') {
        printG = true;
    } else if (answer == 'a') {
        printA = true;
    } else {
        std::cerr << "Did not get it!" << std::endl;
        return 1;
    }
    printf("Perform downfolding (y/n)?");
    if (scanf (" %c",&answer) != 1) {printf("failed to read a value.\n");}
    if (answer == 'y') {
        doDownfolding = true;
        printf("Size of blocks along the diagonal (full inversion within each block):");
        if (scanf("%i",&blocksize) != 1) {printf("failed to read a value.\n");}
    }

    printf("What to print: real part (r), imaginary part (i), or both (b), or the absolute value (a)?");
    if (scanf (" %c",&answer) != 1) {printf("failed to read a value.\n");}
    bool printReal = false;
    bool printImag = false;
    bool printAbs = false;
    if (answer == 'r') {
        printReal = true;
    } else if (answer == 'i') {
        printImag = true;
    } else if (answer == 'b') {
        printReal = true;
        printImag = true;
    } else if (answer == 'a') {
        printReal = false;
        printImag = false;
        printAbs = true;
    } else {
        std::cerr << "Did not get it!" << std::endl;
        return 1;
    }

    std::string line;
    std::vector<std::complex<double> > Hk;
    std::vector<double> kPoints;
    int nkpoints, dim;
    double e_real, e_imag;

//  read k-points and H(k) from a Wien2K *.hk file
    std::ifstream hkfile;
    hkfile.open(infile, std::ifstream::in);
    if (hkfile) {
        getline(hkfile, line);
        nkpoints = atoi (line.substr(1,12).c_str() );
        dim = atoi (line.substr(12,7).c_str() );
        while (getline(hkfile, line)){
       	    for (int x=0; x < 3; ++x) {
                kPoints.push_back(atof (line.substr(2 + 12*x,12).c_str() ) );
            }
	    for (int x=0; x < dim; ++x){
                getline(hkfile, line);
	        for (int y=0; y < dim; ++y){
	      	    e_real = atof(line.substr( 2 + 42*y, 19).c_str());
	      	    e_imag = atof(line.substr(23 + 42*y, 19).c_str());
	       	    Hk.push_back(e_real + e_imag * j);
       	        }
            }
  	}
	hkfile.close();
    } else {
        std::cerr << "Can not read the *.hk file." << std::endl;
        return 1;
    }

//  check whether `dim` is a multiple of `blocksize`
    if ((doDownfolding == true) && (dim % blocksize != 0)) {
        std::cerr << "ERROR: blocksize(" << blocksize
                  << ") must be a multiple of dim(" << dim << ")" << std::endl;
        return 1;
    }

//  calculate the Matsubara frequencies
    std::vector<double> Iw;
    for (int wn=0; wn < nw; ++wn) {Iw.push_back((2*wn + 1) * pi / beta);};
    

//  ACTUAL CALCULATIONS
//  Calculate the Green's function G and the hybridization function F: 
    std::vector<std::complex<double> > G, F;
    for (int windex=0; windex < nw; ++windex){
        float wstep = (wmax - wmin) / (nw - 1);
	std::complex<double> w = (iw==true) ? j*Iw[windex] : wmin + windex*wstep + idelta;
        Eigen::MatrixXcd Gw = Eigen::MatrixXcd::Zero(dim, dim);
	Eigen::MatrixXcd Fw;
        for (unsigned k=0; k < Hk.size()/dim/dim; ++k){
            Eigen::MatrixXcd Gwkinv = Eigen::MatrixXcd::Zero(dim, dim);
            for (int ii = 0; ii < dim; ++ii){
	        for (int jj = 0; jj < dim; ++jj){
	            Gwkinv(ii,jj)  = (ii==jj) ? w + mu : 0; 
                    Gwkinv(ii,jj) -= Hk[k*dim*dim + ii*dim + jj];
	        }
	    }
//            if (doBlockWise){
//                Eigen::MatrixXcd Gw_blockdiag = Eigen::MatrixXcd::Zero(dim, dim);
//                downfold(&Gwkinv, &Gw_blockdiag, blocksize);
//                Gw += Gw_blockdiag;
//            } else {
//	        Gw += (Gwkinv.inverse());
//            }
            if (doDownfolding) {
                Eigen::MatrixXcd Gwk_blockdiag = Eigen::MatrixXcd::Zero(dim, dim);
                Eigen::MatrixXcd Gwk = (Gwkinv.inverse());
                downfold(&Gwk, &Gwk_blockdiag, blocksize);
                Gw += Gwk_blockdiag;
            } else {
                Gw += (Gwkinv.inverse());
            }
        }
        Gw /= Hk.size()/dim/dim;
        Fw = (w + mu) * Eigen::MatrixXcd::Identity(dim,dim) - (Gw.inverse());
        for (int ii=0; ii < dim; ++ii){
            for (int jj=0; jj < dim; ++jj) {
                G.push_back( Gw(ii,jj) );
                F.push_back( Fw(ii,jj) );
            }
        }
    }

//  OUTPUT
//  print the header
    std::ofstream output;
    output.open (outfile, std::ofstream::out | std::ofstream::trunc);

    output << "# Number of k-points: " << nkpoints << std::endl;
    output << "# Hamiltonian matrix: " << dim << "x" << dim << std::endl;
    output << "#" << std::endl << "# G(w) = Sum_k [ w*1 + mu*1 - H(k) ]^-1"
           << std::endl;
    output << "# Delta(w) = w + mu - G(w)^-1   (we assume Sigma(w) = 0 + 0j)"
           << std::endl << "#" << std::endl ;

//  print the legend
    output << "# Column 1: frequency (eV)" << std::endl;
    int colcount = 2; 	
    if (printDiag){
        for (int ii = 0; ii < dim; ++ii){
            if (printDelta){
                if (printReal){
                    output << "# Column " << colcount << ": "
                           << "Delta[" << ii + 1 << "][" << ii + 1
                           << "], real part" << std::endl;
                    ++colcount;
                }
                if (printImag){
                    output << "# Column " << colcount << ": "
                           << "Delta[" << ii + 1 << "][" << ii + 1
                           << "], imaginary part" << std::endl;
                    ++colcount;
                }
                if (printAbs){
                    output << "# Column " << colcount << ": "
                           << "Delta[" << ii + 1 << "][" << ii + 1
                           << "], absolute value" << std::endl;
                    ++colcount;
                }
            }
            if (printG){
                if (printReal){
                    output << "# Column " << colcount << ": "
                           << "G[" << ii + 1 << "][" << ii + 1
                           << "], real part" << std::endl;
                    ++colcount;
                }
                if (printImag){
                    output << "# Column " << colcount << ": "
                           << "G[" << ii + 1 << "][" << ii + 1
                           << "], imaginary part" << std::endl;
                    ++colcount;
                }
                if (printAbs){
                    output << "# Column " << colcount << ": "
                           << "G[" << ii + 1 << "][" << ii + 1
                           << "], absolute value" << std::endl;
                    ++colcount;
                }
            }
            if (printA){
                output << "# Column " << colcount << ": "
                       << "A[" << ii + 1 << "][" << ii + 1
                       << "]" << std::endl;
                ++colcount;
            }
        }
    }

    if (printNonDiag){
        for (int ii=0; ii < dim; ++ii){
            for (int jj=ii+1; jj < dim; ++jj){
                if (printDelta){
                    if (printReal){
                        output << "# Column " << colcount << ": " 
                               << "Delta[" << ii + 1 << "][" << jj + 1
                               << "], real part" << std::endl;
                        ++colcount;
                    }
                    if (printImag){
                        output << "# Column " << colcount << ": "
                               << "Delta[" << ii + 1 << "][" << jj + 1
                               << "], imaginary part" << std::endl;
                        ++colcount;
                    }
                    if (printAbs){
                        output << "# Column " << colcount << ": "
                               << "Delta[" << ii + 1 << "][" << jj + 1
                               << "], absolute value" << std::endl;
                        ++colcount;
                    }
                }
                if (printG){
                    if (printReal){
                        output << "# Column " << colcount << ": "
                               << "G[" << ii + 1 << "][" << jj + 1
                               << "], real part" << std::endl;
                        ++colcount;
                    }
                    if (printImag){
                        output << "# Column " << colcount << ": "
                               << "G[" << ii + 1 << "][" << jj + 1
                               << "], imaginary part" << std::endl;
                        ++colcount;
                    }
                    if (printAbs){
                        output << "# Column " << colcount << ": "
                               << "G[" << ii + 1 << "][" << jj + 1
                               << "], absolute value" << std::endl;
                        ++colcount;
                    }
                }
            }
        }
    }

    for (int windex=0; windex < nw; ++windex){

        float wstep = (wmax - wmin) / (nw - 1);
        std::complex<double> w = wmin + windex * wstep; 
        output << ( (iw==true) ? Iw[windex] : w.real() );

        Eigen::MatrixXcd Gw, Fw;
        Gw = Eigen::MatrixXcd::Map(&G[windex*dim*dim], dim, dim);
        Fw = Eigen::MatrixXcd::Map(&F[windex*dim*dim], dim, dim);

//      print diagonal elements + the spectral function
        if (printDiag){
            for (int ii = 0; ii < dim; ++ii){
                if (printDelta){
                    if (printReal) {output << " " << Fw(ii,ii).real();}
                    if (printImag) {output << " " << Fw(ii,ii).imag();}
                    if (printAbs)  {output << " " << abs(Fw(ii,ii));}
                }
                if (printG){
                    if (printReal) {output << " " << Gw(ii,ii).real();}
                    if (printImag) {output << " " << Gw(ii,ii).imag();}
                    if (printAbs)  {output << " " << abs(Gw(ii,ii));}
                }
                if (printA){
                    output << " " << -1. / pi * Gw(ii,ii).imag();
                }
            }
        }

//      print nondiagonal elements 
        if (printNonDiag){
            for (int ii = 0; ii < dim; ++ii){
                for (int jj = ii + 1; jj < dim; ++jj){
    	            if (printDelta){
                        if (printReal) {output << " " << Fw(ii,jj).real();}
                        if (printImag) {output << " " << Fw(ii,jj).imag();}
                        if (printAbs)  {output << " " << abs(Fw(ii,jj));}
                    }
    		    if (printG){
                        if (printReal) {output << " " << Gw(ii,jj).real();}
                        if (printImag) {output << " " << Gw(ii,jj).imag();}
                        if (printAbs)  {output << " " << abs(Gw(ii,jj));}
                    }
                }
            }
        }
        output << std::endl;
    }
    output.close();

    return 0;
}
