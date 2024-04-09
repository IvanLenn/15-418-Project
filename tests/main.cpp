#include <iostream>
#include <vector>
#include <unistd.h>
#include <string>
#include <fstream>
#include "linear_programming_seq.h"

int main(int argc, char *argv[]) {
    std::string input_filename;
    int opt;
    bool default_input = true;
	while ((opt = getopt(argc, argv, "t:f")) != -1) {
		switch (opt) {
		case 'f':
			input_filename = optarg;
        	default_input = false;
        	break;
    	case 't':
        	break;
		default:
			std::cerr << "Usage: " << argv[0] << " -f input_filename [-t any]\n";
			exit(EXIT_FAILURE);
		}
	}

    if (default_input) {
        std::ifstream fin("../tests/test_data/lp_seq.ans");
        if (!fin) {
			std::cerr << "Unable to open file: " << "../tests/test_data/lp_seq.ans" << ".\n";
			exit(EXIT_FAILURE);
		}
        int N = 0;
        fin >> N;
        std::vector<LinearProgrammingAnswer> Ans;
        for (int i = 0; i < N; i++) {
            int flag;
            double answer = 0.0f;
            fin >> flag >> answer;
            if (flag == 0) {
                LinearProgrammingAnswer tmp;
                Ans.push_back(tmp);
                continue;
            }
            if (flag == 2) {
                LinearProgrammingAnswer tmp(0, std::vector<double>(), LinearProgrammingAnswer::Status::Unbounded);
                Ans.push_back(tmp);
                continue;
            }
            LinearProgrammingAnswer tmp(answer, std::vector<double>(), LinearProgrammingAnswer::Status::Bounded);
            Ans.push_back(tmp);
        }
        fin.close();

        fin.open("../tests/test_data/lp_seq.in");
        if (!fin) {
			std::cerr << "Unable to open file: " << "../tests/test_data/lp_seq.in" << ".\n";
			exit(EXIT_FAILURE);
		}
        fin >> N;
        for (int i = 0; i < N; i++) {
            int numvar, numcons;
            fin >> numvar >> numcons;
            std::vector<std::vector<double>> cons;
            for (int j = 0; j < numcons; j++) {
                std::vector<double> con;
                for (int k = 0; k < numvar + 1; k++) {
                    double tmp;
                    fin >> tmp;
                    con.push_back(tmp);
                }
                cons.push_back(con);
            }
            std::vector<double> target;
            for (int j = 0; j < numvar; j++) {
                double tmp;
                fin >> tmp;
                target.push_back(tmp);
            }
            LinearProgrammingSeq lps(numvar);
            lps.AddTarget(target);
            lps.AddCons(cons);
            LinearProgrammingAnswer Tmp = lps.Solve();
            lps.Check();
            if (!(Tmp == Ans[i])) {
                Tmp.Print();
                Ans[i].Print();
                std::cout << "Test " << i << " failed\n";
            }
            else {
                std::cout << "Test " << i << " passed\n";
            }
            
        }
        fin.close();
        return 0;
    }
}