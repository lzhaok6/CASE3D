#include "stdafx.h"
#include "Header.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

struct LOCAL_NODEstruct LOCAL_NODE(int n) {
	int CN;
	int IC;  //loop counter
	int i, j, k, m;
	LOCAL_NODEstruct t;
	t.LNA = new int**[n + 1];
	for (i = 0; i < n + 1;i++) {
		t.LNA[i] = new int*[n + 1];
		for (j = 0; j < n + 1;j++) {
			t.LNA[i][j] = new int[n + 1];
		}
	}

	for (i = 0; i < n + 1; i++) {
		for (j = 0; j < n + 1; j++) {
			for (k = 0; k < n + 1; k++) {
				t.LNA[i][j][k] = 0;
			}
		}
	}

	if (n == 1) {
		t.LNA[0][0][0] = 1;
		t.LNA[n][0][0] = 2;
		t.LNA[n][n][0] = 3;
		t.LNA[0][n][0] = 4;
		t.LNA[0][0][n] = 5;
		t.LNA[n][0][n] = 6;
		t.LNA[n][n][n] = 7;
		t.LNA[0][n][n] = 8;
	}

	if (internalLNA == 1) {
		if (n > 1) {
			for (i = 1; i < n; i++) { //from the 2nd layer to the one before last layer for the nodes on edge 
				t.LNA[0][0][i] = 8 + 4 * (i - 1) + 1;
				t.LNA[n][0][i] = 8 + 4 * (i - 1) + 2;
				t.LNA[n][n][i] = 8 + 4 * (i - 1) + 3;
				t.LNA[0][n][i] = 8 + 4 * (i - 1) + 4;
			}

			for (i = 0; i < n + 1; i++) { //in the z direction from the first layer to the last layer
				for (IC = 0; IC < n - 1; IC++) { //n-1 internal points on the edge
					t.LNA[IC + 1][0][i] = 4 * (n + 1) + 4 * (n - 1) * i + 4 * IC + 1;
					//4 * (n + 1) is the corner nodes on layers
					//4 * (n - 1) * i is the interior nodes added by previous layers
					t.LNA[n][IC + 1][i] = t.LNA[IC + 1][0][i] + 1;
					t.LNA[(n + 1) - (IC + 1) - 1][n][i] = t.LNA[IC + 1][0][i] + 2;
					t.LNA[0][(n + 1) - (IC + 1) - 1][i] = t.LNA[IC + 1][0][i] + 3;
				}
			}
			//finish the node labeling on boundary of every layer, already ((n+1)*4-4)*(n+1) = 4*n*(n+1) points

			if (n == 2) {           //if n==2, there is a node on center
				t.LNA[1][1][0] = 25;
				t.LNA[1][1][1] = 26;
				t.LNA[1][1][2] = 27;
			}

			else {                  //if n>2
				for (i = 0; i < n + 1; i++) { //label the internal nodes layer by layer
					k = 1;  //interior layer counter
					CN = n - 1;  //the number of node per side of first interior layer

					while (CN >= 2) {
						//std::cout << "the value of CN is:  " << CN << std::endl;
						//std::cout << "the value of k is:  " << k << std::endl;
						if (k == 1) {
							t.LNA[k][k][i] = 4 * n * (n + 1) + (n - 1)*(n - 1)*i + 1;
							//4*n*(n+1) is the points previously labeled
							//(n - 1)*(n - 1)*i is the interior node number from previous layers  
							//when k=1, no additional node from previous interior layer 
							//std::cout << "test initial node:   " << t.LNA[k][k][i] << std::endl;
						}
						else { //starting from k=2
							//std::cout << "the loop is entered" << std::endl;
							//t.LNA[k][k][i] = 4 * n * (n + 1) + (n - 1)*(n - 1)*i + 4 * (n + 3 - 2 * k) + 1;
							//4*(n+3-2*k) is the additional node from previous interior layer, starting from k=2 
							t.LNA[k][k][i] = 1 + t.LNA[k - 1][k][i];
							//std::cout << "test initial node:   "<<t.LNA[k - 1][k][i] << std::endl;
						}
						//everything else would be fine if LNA[k][k][i] is found as the reference point for the current layer
						t.LNA[n - k][k][i] = t.LNA[k][k][i] + 1;
						t.LNA[n - k][n - k][i] = t.LNA[k][k][i] + 2;
						t.LNA[k][n - k][i] = t.LNA[k][k][i] + 3;

						//interior nodes on interior layers
						if (CN > 0) {
							for (j = 0; j < CN - 2; j++) {
								t.LNA[1 + k + (j + 1) - 1][k][i] = t.LNA[k][k][i] + 4 * (j + 1);
								t.LNA[n - k][1 + k + (j + 1) - 1][i] = t.LNA[1 + k + j][k][i] + 1;
								t.LNA[n + 1 - k - (j + 1) - 1][n - k][i] = t.LNA[1 + k + j][k][i] + 2;
								t.LNA[k][n + 1 - k - (j + 1) - 1][i] = t.LNA[1 + k + j][k][i] + 3;
							}
						}
						k = k + 1;
						//std::cout << "the value of k is:  " << k << std::endl;
						CN = CN - 2; //the number of node per side of next interior layer shrink by 2 nodes 

						//check for center node
						if (CN == 1) {
							t.LNA[k][k][i] = 4 * n * (n + 1) + (n - 1)*(n - 1) * i + (n - 1)*(n - 1);
							//(n - 1)*(n - 1) - 1 is the added points other than center node
						}
						//proceed to next interior layer
					}
					//end of while loop
				}
			}
		}
	}
	//import mesh from outside
	else {
		if (n > 1) {
			int u, v, w;
			int NNODE;
			int *ELEMENT_POINT;
			double **FGCOORD;
			std::string filename = "LNA_info_N=" + std::to_string(N) + ".txt";
			std::ifstream infile(filename);   //LNA_info_N=x.txt is for xth order LNA
			if (!infile) {
				std::cerr << "can not open the LNA information file" << std::endl;
				exit(1);
			}
			//read the file 
			int nodestart = 0;
			int nodeend = 0;
			int ele_line = 0;
			std::vector<std::vector<std::string>> output;
			int ct = -1;
			std::string csvLine;
			while (getline(infile, csvLine))
			{
				ct = ct + 1; //the current line number (starting from 0)
				std::istringstream csvStream(csvLine);
				std::vector<std::string> csvColumn;
				std::string csvElement;
				while (getline(csvStream, csvElement, ' '))
				{
					csvColumn.push_back(csvElement);
				}
				output.push_back(csvColumn);

				if (csvColumn[0] == "$Nodes") {
					nodestart = ct + 2; //the node starts from the next line
				}
				if (csvColumn[0] == "$EndNodes") {
					nodeend = ct - 1; //the node starts from the next line
				}
				if (csvColumn[0] == "$Elements") {
					ele_line = ct + 2; //the line corresponding to element connectivity
				}
			}

			NNODE = stoi(output[nodestart - 1][0]);
			ELEMENT_POINT = new int[NINT*NINT*NINT]; //Connectivity matrix
			for (i = 0; i < NINT*NINT*NINT; i++) {
				ELEMENT_POINT[i] = 0;
			}
			for (i = 0; i < NINT*NINT*NINT; i++) {
				ELEMENT_POINT[i] = stoi(output[ele_line][6 + i]);
			}

			FGCOORD = new double*[NNODE];
			for (i = 0; i < NNODE; i++) {
				FGCOORD[i] = new double[3];
			}

			for (i = 0; i < NNODE; i++) {
				for (j = 0; j < 3; j++) {
					FGCOORD[i][j] = stod(output[nodestart + i][1 + j]);
				}
			}

			for (i = 0; i < NINT*NINT*NINT; i++) {
				for (j = 0; j < 3; j++) {
					t.L[i][j] = 0;
				}
			}

			for (i = 0; i < NINT*NINT*NINT; i++) {
				u = round(FGCOORD[ELEMENT_POINT[i] - 1][0] / (1.0 / N));
				v = round(FGCOORD[ELEMENT_POINT[i] - 1][1] / (1.0 / N));
				w = round(FGCOORD[ELEMENT_POINT[i] - 1][2] / (1.0 / N));
				t.L[i][0] = u; t.L[i][1] = v; t.L[i][2] = w;
				t.LNA[u][v][w] = i + 1;
			}
		}
	}
	return t;
}