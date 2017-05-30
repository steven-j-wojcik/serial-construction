#include "graph.h"

int main(int argc, char * argv[]) 
{

	double backstart = MPI_Wtime();
	G.Backtrack(-1);
	double backend = MPI_Wtime() - backstart;
	double constructStart = MPI_Wtime();
	G.Construction(trialsPP); 
	G.Calcval(trialsPP);
	double constructEnd = MPI_Wtime() - constructStart;
	G.setEstimatedRel();
	G.setExactMinCutPath();
	G.setEstimMinCutPath();
	G.Print("output.txt");
	ofstream outFile;
	outFile.open("dodec8np.txt", std::ios::app);
	outFile << "=========================" << endl;
	outFile << endl << "Number of Processors: " << numProcs << endl;
	outFile << "Number of Trials per processor: " << trialsPP << endl;
	outFile << "(proc 0 only does back tracking, doesnt count towards trials)" << endl;
	outFile << "Total number of trials: " << trialsPP << endl;
	outFile << "Backtracking time: " << backend << endl;
	outFile << "Construction time: " << constructEnd << endl;
	outFile << "Constructoin time is the time to aggragte all of the pathsets/totalTrials" << endl;
	outFile.close();

	return 0;
} 