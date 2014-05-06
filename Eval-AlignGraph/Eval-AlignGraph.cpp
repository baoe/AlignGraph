//**********************************************************************************
//* Title: Eval-AlignGraph: Evaluation tool distributed with AlignGraph
//* Platform: 64-Bit Linux
//* Author: Ergude Bao
//* Affliation: Department of Computer Science & Engineering
//* University of California, Riverside
//* Date: 03/24/2011
//* Copy Right: Artistic License 2.0
//**********************************************************************************

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <string>
#include <iomanip>
using namespace std;

#define IDENTITY 0.1
#define CUTOFF 1000
#define SIZE 1000000
#define FASTMAP 0

typedef struct posStruct
{
        int targetID;
        unsigned int sourceStart;
        unsigned int sourceEnd;
        unsigned int targetStart;
        unsigned int targetEnd;
        unsigned int sourceGap;
        unsigned int targetGap;
        int fr;
        int alignedBases;
} pos;

void parseBLAT(string buf, int & targetID, unsigned int & targetStart, unsigned int & targetEnd, unsigned int & targetGap, int & sourceID, unsigned int & sourceStart, unsigned int & sourceEnd, unsigned int & sourceGap, unsigned int & sourceSize, int & fragID, int & fr, int & alignedBases)
{
        char targetIDBuf[100] = {'\0'}, targetStartBuf[100] = {'\0'}, targetEndBuf[100] = {'\0'}, targetGapBuf[100] = {'\0'}, sourceIDBuf[100] = {'\0'}, sourceStartBuf[100] = {'\0'}, sourceEndBuf[100] = {'\0'}, sourceGapBuf[100] = {'\0'}, sourceSizeBuf[100] = {'\0'}, blockBuf[100] = {'\0'}, sourceIDBuf1[100] = {'\0'}, sourceIDBuf2[100] = {'\0'};
        int item = 0, i, j = 0, k, ab = 0;
	vector<int> seg;

        for(i = 0; i < buf.size(); i ++)
        {
                if(buf[i] == '	')
                {
                        item ++;
			j = 0;
                        continue;
                }
                if(item == 13)
                        targetIDBuf[j ++] = buf[i];
                if(item == 15)
                        targetStartBuf[j ++] = buf[i];
                if(item == 16)
                        targetEndBuf[j ++] = buf[i];
                if(item == 7)
                        targetGapBuf[j ++] = buf[i];
                if(item == 9)
                        sourceIDBuf[j ++] = buf[i];
                if(item == 11)
                        sourceStartBuf[j ++] = buf[i];
                if(item == 12)
                        sourceEndBuf[j ++] = buf[i];
                if(item == 5)
                        sourceGapBuf[j ++] = buf[i];
                if(item == 10)
                        sourceSizeBuf[j ++] = buf[i];
		if(item == 8)
			fr = buf[i] == '+' ? 0 : 1;
                if(item == 18)
                {
                        if(buf[i] == ',')
                        {
                                seg.push_back(atoi(blockBuf));

                                for(k = 0; k < j; k ++)
                                        blockBuf[k] = '\0';
                                j = 0;
                        }
                        else
                                blockBuf[j ++] = buf[i];
                }
        }
        targetID = atoi(targetIDBuf); //IMPORTANT! HAS TO BE ADJUSTED FOR MORE THAN ONE CHRS
//	if(targetID > 10) targetID = 0;
        targetStart = atoi(targetStartBuf);
        targetEnd = atoi(targetEndBuf);
        targetGap = atoi(targetGapBuf);
        sourceGap = atoi(sourceGapBuf);
        sourceSize = atoi(sourceSizeBuf);

	for(i = 0; i < 100; i ++)
		if(sourceIDBuf[i] == '.')
			break;
	if(i == 100)
	{
		sourceID = atoi(sourceIDBuf);
		fragID = 0;
	}
	else
	{
		for(j = 0; j < i; j ++)
			sourceIDBuf1[j] = sourceIDBuf[j];
		for(j = i + 1; j < 100; j ++)
			sourceIDBuf2[j - i - 1] = sourceIDBuf[j];
		sourceID = atoi(sourceIDBuf1);
		fragID = atoi(sourceIDBuf2);
	}
	sourceStart = fragID * SIZE + atoi(sourceStartBuf);
	sourceEnd = fragID * SIZE + atoi(sourceEndBuf);

	for(i = 0; i < seg.size(); i ++)
		ab = ab + seg[i];
	alignedBases = ab;
}

int conflict(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)
{
        if(x1 <= x2 && x2 <= y1 && y1 <= y2 && (int)y1 - (int)x2 >= 100 || x2 <= x1 && x1 <= y2 && y2 <= y1 && (int)y2 - (int)x1 >= 100 || x1 <= x2 && x2 <= y2 && y2 <= y1 && (int)y2 - (int)x2 >= 100 || x2 <= x1 && x1 <= y1 && y1 <= y2 && (int)y1 - (int)x1 >= 100 ||
	x1 <= x2 && y2 <= y1 || x2 <= x1 && y1 <= y2)
                return 1;
        else
                return 0;
}

int close(unsigned int y1, unsigned int x2, unsigned int threshold)
{
	if(abs((int)x2 - (int)y1) < threshold)
		return 1;
	else
		return 0;
}

vector<vector<int> > loadGenome()
{
	ifstream in;
	string buf;
	vector<int> gb;
	int i;
	vector<vector<int> > genome;

	in.open("etmp/_genome.fa");
	if(in.is_open())
	{
		while(in.good())
		{
			getline(in, buf);
			if(buf[0] == 0) break;

			if(buf[0] == '>')
				genome.push_back(gb);
			else
				for(i = 0; i < buf.size(); i ++)
					genome[genome.size() - 1].push_back(0);
		}
	}
	else
	{
		cout << "CANNOT OPEN FILE!" << endl;
		exit(-1);
	}

	return genome;
}

vector<vector<char> > loadContigs()
{
	ifstream in;
	string buf;
	int i, j, initIDBak = -1;
	char initIDBuf[100];
	vector<vector<char> > initContigs;
	vector<char> ic;

	in.open("etmp/_contigs.fa");
	if(in.is_open())
	{
		while(in.good())
		{
			getline(in, buf);
			if(buf[0] == 0) break;

			if(buf[0] == '>')
			{
				for(i = 1, j = 0; i < buf.size() && buf[i] != '.'; i ++, j ++) initIDBuf[j] = buf[i];
				for(; j < 100; j ++) initIDBuf[j] = '\0';
				if(atoi(initIDBuf) > initIDBak)
				{
					initContigs.push_back(ic);
					initIDBak = atoi(initIDBuf);
				}
				continue;
			}
			else
				for(i = 0; i < buf.size(); i ++)
					initContigs[initContigs.size() - 1].push_back(buf[i]);
		}
	}
	else
	{
		cout << "CANNOT OPEN FILE!" << endl;
		exit(-1);
	}

	return initContigs;
}

vector<vector<pos> > loadContigsAlignment(int size)
{
	ifstream in;
	vector<vector<pos> > positions;
	vector<pos> position;
	pos p, p0;
	string buf;
	unsigned int targetStart, targetEnd, targetGap, sourceStart, sourceEnd, sourceGap, sourceSize;
	int targetID, sourceID, keep, i, j, k, fragID, fr, alignedBases;

	for(i = 0; i < size; i ++)
		positions.push_back(position);
	p0.targetID = p0.sourceStart = p0.sourceEnd = p0.targetStart = p0.targetEnd = p0.fr = p0.sourceGap = p0.targetGap = p0.alignedBases = -1;

	in.open("etmp/_contigs_genome.psl");
	if(in.is_open())
	{
		while(in.good())
		{
			getline(in, buf);
			if(buf[0] == 0) break;

			parseBLAT(buf, targetID, targetStart, targetEnd, targetGap, sourceID, sourceStart, sourceEnd, sourceGap, sourceSize, fragID, fr, alignedBases);
			if(sourceEnd - sourceStart >= 100 && (double)(sourceEnd - sourceStart - sourceGap) / (sourceEnd - sourceStart) >= IDENTITY && (double)(targetEnd - targetStart - targetGap) / (double)(targetEnd - targetStart) >= IDENTITY)
			{
				keep = 1;
				for(i = 0; i < positions[sourceID].size(); i ++)
					if(positions[sourceID][i].targetID != -1 && targetID == positions[sourceID][i].targetID && conflict(sourceStart, sourceEnd, positions[sourceID][i].sourceStart, positions[sourceID][i].sourceEnd))
					{
						if(sourceEnd - sourceStart < positions[sourceID][i].sourceEnd - positions[sourceID][i].sourceStart)
							keep = 0;
						else
							positions[sourceID][i] = p0;
					}
				if(keep)
				{
					p.sourceStart = sourceStart;
					p.sourceEnd = sourceEnd;
					p.targetStart = targetStart;
					p.targetEnd = targetEnd;
					p.targetID = targetID;
					p.sourceGap = sourceGap;
					p.targetGap = targetGap;
					p.fr = fr;
					p.alignedBases = alignedBases;
					positions[sourceID].push_back(p);
				}
			}
		}
	}
	else
	{
		cout << "CANNOT OPEN FILE!" << endl;
		exit(-1);
	}

	for(i = 0; i < positions.size(); i ++)
	{
		for(j = 0; j < positions[i].size(); j ++)
		{
			for(k = 0; k < positions[i].size(); k ++)
			{
				if(k != j && positions[i][j].targetID != -1 && positions[i][k].targetID != -1 && positions[i][j].targetID == positions[i][k].targetID && close(positions[i][j].sourceEnd, positions[i][k].sourceStart, abs((int)positions[i][j].sourceEnd - (int)positions[i][j].sourceStart) / 10) && close(positions[i][j].targetEnd, positions[i][k].targetStart, abs((int)positions[i][j].targetEnd - (int)positions[i][j].targetStart) / 10) && positions[i][j].fr == positions[i][k].fr)
				{

					positions[i][j].sourceEnd = positions[i][k].sourceEnd;
					positions[i][j].targetEnd = positions[i][k].targetEnd;
					positions[i][j].sourceGap = positions[i][j].sourceGap + positions[i][k].sourceGap;
					positions[i][j].targetGap = positions[i][j].targetGap + positions[i][k].targetGap;
					positions[i][j].alignedBases = positions[i][j].alignedBases + positions[i][k].alignedBases;
					positions[i][k] = p0;
					k = 0;
				}
			}
		}
	}

	for(i = 0; i < positions.size(); i ++)
		for(j = 0; j < positions[i].size(); j ++)
			for(k = j + 1; k < positions[i].size(); k ++)
			{
				if(positions[i][j].targetID != -1 && positions[i][k].targetID != -1 && conflict(positions[i][j].sourceStart, positions[i][j].sourceEnd, positions[i][k].sourceStart, positions[i][k].sourceEnd) == 1)
				{
					if(positions[i][j].sourceEnd - positions[i][j].sourceStart > positions[i][k].sourceEnd - positions[i][k].sourceStart)
						positions[i][k] = p0;
					else
					{
						positions[i][j] = p0;
						break;
					}
				}
			}
//remove duplicated alignments to different chrs

	return positions;
}

void analyze(char argv[], vector<vector<int> > & genome, vector<vector<char> > & initContigs,  vector<vector<pos> > & positions)
{
	ofstream out;
	int max, misassembly, i, j, k, alignedBases, totalBases, errors, sum, coveredLength, contigBases, totalLength;
	vector<int> trueContigLengths;
	vector<double> identity;
	double totalIdentity;

	max = misassembly = 0;
	for(i = 0; i < positions.size(); i ++)
	{
		for(j = 0; j < positions[i].size(); j ++)
			if(positions[i][j].targetID != -1 && (double)(positions[i][j].sourceEnd - positions[i][j].sourceStart) / initContigs[i].size() >= 0.8)
			{
				trueContigLengths.push_back(positions[i][j].sourceEnd - positions[i][j].sourceStart);

				for(k = positions[i][j].targetStart; k < positions[i][j].targetEnd; k ++)
					genome[positions[i][j].targetID][k] = 1;

				if(max < positions[i][j].sourceEnd - positions[i][j].sourceStart)
					max = positions[i][j].sourceEnd - positions[i][j].sourceStart;

				alignedBases = positions[i][j].alignedBases;
				totalBases = positions[i][j].targetEnd - positions[i][j].targetStart + positions[i][j].targetGap;
				identity.push_back((double) alignedBases * trueContigLengths[trueContigLengths.size() - 1] / totalBases);

				goto end;
			}

		for(j = 0, errors = 0; j < positions[i].size(); j ++)
		{
			if(positions[i][j].targetID != -1)
			{
				trueContigLengths.push_back(positions[i][j].sourceEnd - positions[i][j].sourceStart);

				for(k = positions[i][j].targetStart; k < positions[i][j].targetEnd; k ++)
					genome[positions[i][j].targetID][k] = 1;

				if(max < positions[i][j].sourceEnd - positions[i][j].sourceStart)
					max = positions[i][j].sourceEnd - positions[i][j].sourceStart;

				alignedBases = positions[i][j].alignedBases;
				totalBases = positions[i][j].targetEnd - positions[i][j].targetStart + positions[i][j].targetGap;
				identity.push_back((double) alignedBases * trueContigLengths[trueContigLengths.size() - 1] / totalBases);

				errors ++;
			}
		}
		if(errors == 0)
			;//misassembly ++;//do nothing
		else if(errors == 1)
			misassembly ++;
		else if(errors >= 2)
			misassembly = misassembly + errors - 1;
end:;
	}

	out.open(argv);
	if(out.is_open() == 0) {cout << "CANNOT OPEN FILE!"; exit(-1);}
	out << setw(21) << left << "#contigs" << initContigs.size() << endl;
	out << setw(21) << left << "#true contigs" << trueContigLengths.size() << endl;

	for(i = 0, totalLength = 0; i < trueContigLengths.size(); i ++)
		totalLength = totalLength + trueContigLengths[i];
	sort(trueContigLengths.begin(), trueContigLengths.end());
	for(i = trueContigLengths.size() - 1, sum = 0; i >= 0; i --)
	{
		sum = sum + trueContigLengths[i];
		if(sum > totalLength / 2) break;
	}
	out << setw(21) << left << "N50" << trueContigLengths[i] << endl;

	for(i = 0, coveredLength = 0; i < genome.size(); i ++)
		for(j = 0; j < genome[i].size(); j ++)
			if(genome[i][j]) coveredLength ++;
	out << setw(21) << left << "covered length" << coveredLength << endl;

	for(i = 0, contigBases = 0; i < initContigs.size(); i ++)
		contigBases = contigBases + initContigs[i].size();

	out << setw(21) << left << "average length" << totalLength / trueContigLengths.size() << endl;

	out << setw(21) << left << "maximum length" << max << endl;

	out << setw(21) << left << "MPMB" << (double) misassembly / ((double) (contigBases) / 1000000) << endl;

	for(i = 0, totalIdentity = 0; i < identity.size(); i ++)
		totalIdentity = totalIdentity + identity[i];
	out << setw(21) << left << "average identity" << totalIdentity / totalLength << endl;
}

void formalizeGenome(char argv[])
{
	ifstream in;
	ofstream out;
	string buf;
	vector<vector<char> > genome;
	vector<char> g;
	int i, j;

	in.open(argv);
	if(in.is_open())
	{
		while(in.good())
		{
			getline(in, buf);
			if(buf[0] == 0) break;

			if(buf[0] == '>')
				genome.push_back(g);
			else
				for(i = 0; i < buf.size(); i ++)
					genome[genome.size() - 1].push_back(buf[i]);
		}
	}
	else
	{
		cout << "CANNOT OPEN FILE!" << endl;
		exit(-1);
	}

	out.open("etmp/_genome.fa");
	if(out.is_open())
	{
		for(i = 0; i < genome.size(); i ++)
		{
			out << ">" << i << endl;
			for(j = 0; j < genome[i].size(); j ++)
			{
				out << genome[i][j];
				if(j == genome[i].size() - 1 || (j + 1) % 60 == 0)
					out << endl;
			}
		}
	}
	else
	{
		cout << "CANNOT OPEN FILE!" << endl;
		exit(-1);
	}
}

void formalizeContigs(char argv[])
{
	ifstream in;
	ofstream out;
	string buf;
	vector<vector<char> > contigs;
	vector<char> c;
	vector<vector<char> > contigsBuf;
	int i, j, k, seqID = 0;

	in.open(argv);
	if(in.is_open())
	{
		while(in.good())
		{
			getline(in, buf);
			if(buf[0] == 0) break;

			if(buf[0] == '>')
				contigs.push_back(c);
			else
				for(i = 0; i < buf.size(); i ++)
					contigs[contigs.size() - 1].push_back(buf[i]);
		}
	}
	else
	{
		cout << "CANNOT OPEN FILE!" << endl;
		exit(-1);
	}

	out.open("etmp/_contigs.fa");
	if(out.is_open())
	{
		for(i = 0; i < contigs.size(); i ++)
		{
			if(contigs[i].size() >= CUTOFF)
			{
				if(contigs[i].size() < SIZE)
				{
					out << ">" << seqID ++ << endl;
					for(j = 0; j < contigs[i].size(); j ++)
					{
						out << contigs[i][j];
						if((j + 1) % 60 == 0 || j == contigs[i].size() - 1)
							out << endl;
					}
				}
				else
				{
					contigsBuf.clear();
					contigsBuf.push_back(c);
					for(j = 0; j < contigs[i].size(); j ++)
					{
						contigsBuf[contigsBuf.size() - 1].push_back(contigs[i][j]);
						if((j + 1) % SIZE == 0 && j < contigs[i].size() - 1) contigsBuf.push_back(c);
					}
					
					for(j = 0; j < contigsBuf.size(); j ++)
					{
						out << ">" << seqID << "." << j << endl;
						for(k = 0; k < contigsBuf[j].size(); k ++)
						{
							out << contigsBuf[j][k];
							if((k + 1) % 60 == 0 || k == contigsBuf[j].size() - 1)
								out << endl;
						}
					}

					seqID ++;
				}
			}
		}
	}
	else
	{
		cout << "CANNOT OPEN FILE!" << endl;
		exit(-1);
	}
}

void makeAlignment()
{
	if(FASTMAP == 1)
		system("blat etmp/_genome.fa etmp/_contigs.fa -noHead etmp/_contigs_genome.psl -fastMap > blat_doc.txt");
	else
		system("blat etmp/_genome.fa etmp/_contigs.fa -noHead etmp/_contigs_genome.psl > blat_doc.txt");
}

void print()
{
	cout << "Eval-AlignGraph arg1, arg2, arg3" << endl;
	cout << "arg1 = file of target genome in FASTA format" << endl;
	cout << "arg2 = file of contigs in FASTA format" << endl;
	cout << "arg3 = file name for outputted statistics" << endl;
}

int main(int argc, char * argv[])
{
        vector<vector<pos> > positions;
        vector<vector<int> > genome;
        vector<vector<char> > initContigs;

	cout << "Eval-AlignGraph: Evaluation tool distributed with AlignGraph" << endl;
	cout << "By Ergude Bao, CS Department, UC-Riverside. All Rights Reserved" << endl << endl;

        if(argc != 4) {print(); return 0;}

	system("test -d \"etmp\"; t=$?; if [ $t -eq 1 ]; then mkdir etmp; fi");
	formalizeGenome(argv[1]);
	formalizeContigs(argv[2]);
	makeAlignment();
	cout << "(1) Alignment finished" << endl;

        genome = loadGenome();
        initContigs = loadContigs();
        positions = loadContigsAlignment(initContigs.size());
        analyze(argv[3], genome, initContigs, positions);
	cout << "(2) Statistics generated and done :-)" << endl;
}

