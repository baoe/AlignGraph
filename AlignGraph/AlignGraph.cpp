#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <string>
#include <time.h>
using namespace std;

#define TEST
#define OPTIMIZATION

#define MAX 99999
#define INIT_CONTIG_THRESHOLD 0.5//0
#define CONTIG_THRESHOLD 0.5//0.1
#define THRESHOLD 0.6//0.2
#define SD 0//5
#define SI 0//5
#define BATCH 1000000// BATCH = 2n
#define COVERAGE 10
#define EP 5// 5, 20, 100, max
#define LARGE_CHUNK 1000000//20000
#define SMALL_CHUNK 20000

typedef struct structSegment
{
        unsigned int sourceStart;
        unsigned int targetStart;
        unsigned int size;
} Segment;

typedef struct contiMerStruct
{
	char nucleotide;
        unsigned int contigID;
        unsigned int contigOffset;
        //unsigned int previousID;
        //unsigned int previousOffset;
        //unsigned int previousItem;
        unsigned int nextID;
        unsigned int nextOffset;
        unsigned int nextItem;
} ContiMer;

typedef struct nextStruct
{
        unsigned int nextID;
        unsigned int nextOffset;
        unsigned int nextItem;
} Next;

typedef struct previousStruct
{
        unsigned int previousID;
        unsigned int previousOffset;
        unsigned int previousItem;
} Previous;

typedef struct kMerStruct
{
        unsigned int traversed;
        vector<char> s;
        unsigned int contigID;// omitted, since it is rare that two k-mers at the same genome position with the same contig offset have different contig IDs
        unsigned int contigOffset;
        unsigned int contigID0;// omitted with the same reason with above
        unsigned int contigOffset0;
        unsigned int chromosomeID0;
        unsigned int chromosomeOffset0;
        vector<Next> next;
        //vector<Previous> previous;
//      unsigned int nextID;
//      unsigned int nextOffset;
//      unsigned int nextItem;
//      unsigned int previousID;
//      unsigned int previousOffset;
//      unsigned int previousItem;
        int coverage;
	int A, C, G, T, N;
} KMer;

typedef struct baseStruct
{
        vector<KMer> kMer;
        char nucleotide;
        vector<ContiMer> contiMer;
} Base;

vector<vector<Base> > genome;

string itoa(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

void parseBOWTIE(string buf, unsigned int & targetID, unsigned int & targetStart, unsigned int & targetEnd, unsigned int & targetGap, unsigned int & sourceID, unsigned int & sourceStart, unsigned int & sourceEnd, unsigned int & sourceGap, unsigned int & sourceSize, vector<Segment> & seg, unsigned int & fr)
{
//void parseBOWTIE(char buf[], int & targetID, unsigned int & targetStart, unsigned int & targetInsertion, unsigned int & targetDeletion, unsigned int & sourceSize)
//{
	int item = 0, i, j0 = 0, j1 = 0, j2 = 0, j3 = 0, j4 = 0, k, insertion = 0, deletion = 0, total = 0, start = 0, end = 0, tag = 1;
	char sourceIDBuf[1000] = {'\n'}, targetIDBuf[1000] = {'\n'}, sourceStartBuf[1000] = {'\n'}, targetStartBuf[1000] = {'\n'}, CIGARBuf[1000] = {'\n'}, frBuf[1000] = {'\n'};
	Segment s;
        
	for(i = 0; i < buf.size(); i ++)
        {
                if(buf[i] == '	')
                {
                        item ++;
                        continue;
                }

		if(item == 0)
			sourceIDBuf[j3 ++] = buf[i];
		if(item == 1)
			frBuf[j4 ++] = buf[i];
                if(item == 2)
		{
			if(buf[i] == '*')
			{
				sourceID = atoi(sourceIDBuf);
				fr = ((atoi(frBuf) & 0x00000010) == 0x00000010) ? 1 : 0;
				targetID = targetStart = targetEnd = targetGap = sourceStart = sourceEnd = sourceGap = sourceSize = -1;
				return;
			}	
			targetIDBuf[j0 ++] = buf[i];
		}
		if(item == 3)
			targetStartBuf[j1 ++] = buf[i];
		if(item == 5)
			if(buf[i] == '0' || buf[i] == '1' || buf[i] == '2' || buf[i] == '3' || buf[i] == '4' || buf[i] == '5' || buf[i] == '6' || buf[i] == '7' || buf[i] == '8' || buf[i] == '9')
				CIGARBuf[j2 ++] = buf[i];
			else if(buf[i] == 'I')
			{
				insertion = insertion + atoi(CIGARBuf);
				total = total + atoi(CIGARBuf);
				for(k = 0; k < 1000; k ++)
					CIGARBuf[k] = '\n';
				j2 = 0;
			}
			else if(buf[i] == 'D')
			{
				deletion = deletion + atoi(CIGARBuf);
//				total = total + atoi(CIGARBuf);
				for(k = 0; k < 1000; k ++)
					CIGARBuf[k] = '\n';
				j2 = 0;
			}
			else if(buf[i] == 'M')
			{
				s.sourceStart = total;
				s.targetStart = atoi(targetStartBuf) + total + deletion - start - insertion - 1;// fixed a bug; fixed another bug: offset is 1-based
				s.size = atoi(CIGARBuf);
				seg.push_back(s);
				total = total + atoi(CIGARBuf);
				for(k = 0; k < 1000; k ++)
					CIGARBuf[k] = '\n';
				j2 = tag = 0;
			}
			else if(buf[i] == 'S' && tag)
			{
				start = atoi(CIGARBuf);
				total = total + atoi(CIGARBuf);
				for(k = 0; k < 1000; k ++)
					CIGARBuf[k] = '\n';
				j2 = tag = 0;
			}
			else if(buf[i] == 'S')
			{
				end = atoi(CIGARBuf);
				total = total + atoi(CIGARBuf);
				for(k = 0; k < 1000; k ++)
					CIGARBuf[k] = '\n';
				j2 = 0;
			}
			else
			{
				if(buf[i] != '*')
				{
					cout << "unknown character: " << buf[i] << endl;
					exit(-1);
				}
			}
        }
	sourceID = atoi(sourceIDBuf);
	sourceStart = start;
	sourceEnd = total - end;
	sourceGap = insertion;
	sourceSize = total;
	targetID = targetIDBuf[0] == '*' ? -1 : atoi(targetIDBuf);
	targetStart = atoi(targetStartBuf) - 1;//fixed a bug: offset is 1-based
	targetEnd = targetStart + total + deletion;
	targetGap = deletion;
	fr = ((atoi(frBuf) & 0x00000010) == 0x00000010) ? 1 : 0;
}

/*
typedef struct contiMerStruct
{
	unsigned int contigID;
	unsigned int contigOffset;
	//unsigned int previousID;
	//unsigned int previousOffset;
	//unsigned int previousItem;
	unsigned int nextID;
	unsigned int nextOffset;
	unsigned int nextItem;
} ContiMer;

typedef struct nextStruct
{
	unsigned int nextID;
	unsigned int nextOffset;
	unsigned int nextItem;
} Next;

typedef struct previousStruct
{
	unsigned int previousID;
	unsigned int previousOffset;
	unsigned int previousItem;
} Previous;

typedef struct kMerStruct
{
	unsigned int traversed;
	vector<char> s;
        unsigned int contigID;// omitted, since it is rare that two k-mers at the same genome position with the same contig offset have different contig IDs
        unsigned int contigOffset;
        unsigned int contigID0;// omitted with the same reason with above
        unsigned int contigOffset0;
        unsigned int chromosomeID0;
        unsigned int chromosomeOffset0;
	vector<Next> next;
	//vector<Previous> previous;
//	unsigned int nextID;
//	unsigned int nextOffset;
//	unsigned int nextItem;
//	unsigned int previousID;
//	unsigned int previousOffset;
//	unsigned int previousItem;
	int coverage;
} KMer;

typedef struct baseStruct
{
	vector<KMer> kMer;
	char nucleotide;
	vector<ContiMer> contiMer;
} Base;

vector<vector<Base> > genome;
*/

void loadGenome(vector<vector<Base> > & genome, int chromosomeID)
{
	vector<Base> chromosome;
	string buf;
	Base b;
	int gp, cp, i;
	ifstream g;
	string s;

	s = "tmp/_genome." + itoa(chromosomeID) + ".fa";
	g.open(s.c_str());
	if(g.is_open())
	{
		while(g.good())
		{
			getline(g, buf);
			if(buf[0] == 0) break;

			if(buf[0] == '>')
				genome.push_back(chromosome);
			else
				for(i = 0; i < buf.size(); i ++)
				{
					b.nucleotide = buf[i];
					genome[genome.size() - 1].push_back(b);
				}
		}
	}
	else
	{
		cout << "CANNOT OPEN FILE!" << endl;
		exit(-1);
	}
}

typedef struct positionStruct
{
	unsigned int chromosomeID;
	unsigned int chromosomeOffset;
} Position;

typedef struct seqStruct
{
	int adjusted;
	vector<char> nucleotides;
	vector<vector<Position> > positionSets;
	vector<int> frSets;
	int outputted;
	int ID;
} Seq;

void loadSeq(ifstream & in, vector<Seq> & seqs)
{
        string buf;
        int i, sp, ssp;
        Seq s;

        s.adjusted = 0;
	s.outputted = 0;	

        if(in.is_open())
        {
                while(in.good())
                {
                        getline(in, buf);
                        if(buf[0] == 0)
                                break;

                        if(buf[0] == '>')
			{
				for(i = 0; i < buf.size(); i ++)
					if(buf[i] == '.')
						break;
				s.ID = atoi(buf.substr(i + 1, buf.size()).c_str());
                                seqs.push_back(s);
			}
                        else
                        {
                                for(i = 0; i < buf.size(); i ++)
                                        seqs[seqs.size() - 1].nucleotides.push_back(buf[i]);
                        }
                }
        }
        else
        {
                cout << "CANNOT OPEN FILE!" << endl;
                exit(-1);
        }
}

int loadSeq(ifstream & in, vector<Seq> & seqs, int & aliStartID, int & seqStartID)
{
	string buf;
	int i, sp, ssp, p;
	Seq s;

	aliStartID = seqStartID + 1;
	s.adjusted = 0;
	p = 0;
	if(in.is_open())
	{
		while(in.good())
		{
			getline(in, buf);
			if(buf[0] == 0)
				break;

			if(buf[0] == '>')
			{
				seqs.push_back(s);
			}
			else
			{
				for(i = 0; i < buf.size(); i ++)
					seqs[seqs.size() - 1].nucleotides.push_back(buf[i]);

				if(p == 1)
				{
					seqStartID ++;
					if((seqStartID + 1) % BATCH == 0) return 0;
					p = 0;
				}
				else
					p = 1;
			}
		}
		return 1;
	}
	else
	{
		cout << "CANNOT OPEN FILE!" << endl;
		exit(-1);
	}

//        for(sp = 0; sp < seqs.size(); sp ++)
//                for(ssp = 0; ssp < seqs[sp].nucleotides.size(); ssp ++)
//                {
//                        cout << seqs[sp].nucleotides[ssp];
//                        if((ssp + 1) % 60 == 0 || ssp == seqs[sp].nucleotides.size() - 1)
//                                cout << endl;
//                }

}

int parseBLAT(string buf, unsigned int & targetID, unsigned int & targetStart, unsigned int & targetEnd, unsigned int & targetGap, unsigned int & sourceID, unsigned int & sourceStart, unsigned int & sourceEnd, unsigned int & sourceGap, unsigned int & sourceSize, vector<Segment> & seg, unsigned int & fr, unsigned int & targetSize)
{
        char targetIDBuf[100] = {'\0'}, targetStartBuf[100] = {'\0'}, targetEndBuf[100] = {'\0'}, targetGapBuf[100] = {'\0'}, sourceIDBuf[100] = {'\0'}, sourceStartBuf[100] = {'\0'}, sourceEndBuf[100] = {'\0'}, sourceGapBuf[100] = {'\0'}, sourceSizeBuf[100] = {'\0'}, blockBuf[100] = {'\0'}, targetSizeBuf[100] = {'\0'}, realSourceIDBuf[100] = {'\0'}, realBuf[100] = {'\0'};
        int item = 0, i, j0 = 0, j1 = 0, j2 = 0, j3 = 0, j4 = 0, j5 = 0, j6 = 0, j7 = 0, j8 = 0, j9 = 0, j, k, sp, j10 = 0, j11 = 0, j12 = 0;
        Segment s;
	fr = -1;

        seg.clear();

        for(i = 0; i < buf.size(); i ++)
        {
                if(buf[i] == '	')
                {
                        item ++;
			sp = 0;
                        continue;
                }

                if(item == 13)
                        targetIDBuf[j0 ++] = buf[i];
                if(item == 15)
                        targetStartBuf[j1 ++] = buf[i];
                if(item == 16)
                        targetEndBuf[j2 ++] = buf[i];
                if(item == 7)
                        targetGapBuf[j3 ++] = buf[i];
                if(item == 9)
                        sourceIDBuf[j4 ++] = buf[i];
                if(item == 11)
                        sourceStartBuf[j5 ++] = buf[i];
                if(item == 12)
                        sourceEndBuf[j6 ++] = buf[i];
                if(item == 5)
                        sourceGapBuf[j7 ++] = buf[i];
                if(item == 10)
                        sourceSizeBuf[j8 ++] = buf[i];
                if(item == 18)
                {
                        if(buf[i] == ',')
                        {
                                s.sourceStart = s.targetStart = -1;
                                s.size = atoi(blockBuf);
                                seg.push_back(s);

                                for(k = 0; k < j9; k ++)
                                        blockBuf[k] = '\0';
                                j9 = 0;
                        }
                        else
                                blockBuf[j9 ++] = buf[i];
                }
                if(item == 19)
                {
                        if(buf[i] == ',')
                        {
                                seg[sp ++].sourceStart = atoi(blockBuf);

                                for(k = 0; k < j9; k ++)
                                        blockBuf[k] = '\0';
                                j9 = 0;
                        }
                        else
                                blockBuf[j9 ++] = buf[i];
                }
                if(item == 20)
                {
                        if(buf[i] == '\0')
                                break;
                        if(buf[i] == ',')
                        {
                                seg[sp ++].targetStart = atoi(blockBuf);

                                for(k = 0; k < j9; k ++)
                                        blockBuf[k] = '\0';
                                j9 = 0;
                        }
                        else
                                blockBuf[j9 ++] = buf[i];
                }
                if(item == 8 && fr == -1)
                        fr = buf[i] == '+' ? 0 : 1;
		if(item == 14)
                        targetSizeBuf[j10 ++] = buf[i];
        }
        targetID = atoi(targetIDBuf);
        targetStart = atoi(targetStartBuf);
        targetEnd = atoi(targetEndBuf);
        targetGap = atoi(targetGapBuf);
//      sourceID = atoi(sourceIDBuf);
        sourceStart = atoi(sourceStartBuf);
        sourceEnd = atoi(sourceEndBuf);
        sourceGap = atoi(sourceGapBuf);
//      sourceSize = atoi(sourceSizeBuf);
	targetSize = atoi(targetSizeBuf);

	for(i = 0; i < 100; i ++)
		if(sourceIDBuf[i] == '.')
			break;
	if(i < 100)
	{
		for(j = 0; j < i; j ++)
			realSourceIDBuf[j11 ++] = sourceIDBuf[j];
		for(j = i + 1; j < 100; j ++)
			realBuf[j12 ++] = sourceIDBuf[j];
		sourceID = atoi(realSourceIDBuf);
		sourceSize = atoi(sourceSizeBuf) ;
		return atoi(realBuf);
	}
	else
	{
		sourceID = atoi(sourceIDBuf);
		sourceSize = atoi(sourceSizeBuf);
		return sourceSize;
	}	
}

int keepPositions(vector<Seq> & contigs, unsigned int sourceID, vector<Segment> & segs, double threshold)
//This is only used to complete the combination of BLAT's local alignments
{
	int pp, match;

	if(sourceID == -1) return 1;
	if(contigs[sourceID].positionSets.size() == 0) return 1;
	match = 0;
	for(pp = 0; pp < contigs[sourceID].positionSets[contigs[sourceID].positionSets.size() - 1].size(); pp ++)
	{
		if(contigs[sourceID].positionSets[contigs[sourceID].positionSets.size() - 1][pp].chromosomeID != -1)
			match ++;
	}
	if((double) match / contigs[sourceID].positionSets[contigs[sourceID].positionSets.size() - 1].size() >= threshold)// 0.8 to 0 for test
		return 1;
	else
		return 0;
}

void complement(vector<Seq> & contigs, unsigned int sourceID)
{
        unsigned int np;

        for(np = 0; np < contigs[sourceID].nucleotides.size(); np ++)
                if(contigs[sourceID].nucleotides[np] == 'A') contigs[sourceID].nucleotides[np] = 'T';
                else if(contigs[sourceID].nucleotides[np] == 'C') contigs[sourceID].nucleotides[np] = 'G';
                else if(contigs[sourceID].nucleotides[np] == 'G') contigs[sourceID].nucleotides[np] = 'C';
                else if(contigs[sourceID].nucleotides[np] == 'T') contigs[sourceID].nucleotides[np] = 'A';
                //do nothing for 'N'
}

static unsigned int sourceIDBak = -1;
void updateContig(vector<Seq> & contigs, unsigned int sourceID, unsigned int targetID, vector<Segment> & segs, unsigned int fr, double threshold)
{
//	static unsigned int sourceIDBak = -1;
	int sp, ssp, np;
	vector<Position> positions;
	Position p;

	if(targetID == -1) return;
	p.chromosomeID = p.chromosomeOffset = -1;

	if(sourceID != sourceIDBak)
	{
		if(keepPositions(contigs, sourceIDBak, segs, threshold) == 0)
		{
			contigs[sourceIDBak].positionSets.erase(contigs[sourceIDBak].positionSets.end() - 1);
			contigs[sourceIDBak].frSets.erase(contigs[sourceIDBak].frSets.end() - 1);
		}
		contigs[sourceID].positionSets.push_back(positions);
		contigs[sourceID].frSets.push_back(fr);
		for(np = 0; np < contigs[sourceID].nucleotides.size(); np ++)
			contigs[sourceID].positionSets[contigs[sourceID].positionSets.size() - 1].push_back(p);
		//contigs[sourceID].positionSets[contigs[sourceID].positionSets.size() - 1].resize(contigs[sourceID].nucleotides.size(), p);
		sourceIDBak = sourceID;
	}
	else
	{
		for(sp = 0; sp < segs.size(); sp ++)
			for(ssp = segs[sp].sourceStart; ssp < segs[sp].sourceStart + segs[sp].size; ssp ++)
			{
				if(contigs[sourceID].positionSets[contigs[sourceID].positionSets.size() - 1][ssp].chromosomeID != -1)
				{
					if(keepPositions(contigs, sourceID, segs, threshold) == 0)
					{
						contigs[sourceID].positionSets.erase(contigs[sourceID].positionSets.end() - 1);
						contigs[sourceID].frSets.erase(contigs[sourceID].frSets.end() - 1);
					}
					contigs[sourceID].positionSets.push_back(positions);
					contigs[sourceID].frSets.push_back(fr);
					for(np = 0; np < contigs[sourceID].nucleotides.size(); np ++)
						contigs[sourceID].positionSets[contigs[sourceID].positionSets.size() - 1].push_back(p);
					//contigs[sourceID].positionSets[contigs[sourceID].positionSets.size() - 1].resize(contigs[sourceID].nucleotides.size(), p);
					goto cont;
				}
			}
	}

cont:
	for(sp = 0; sp < segs.size(); sp ++)
		for(ssp = 0; ssp < segs[sp].size; ssp ++)
		{
			contigs[sourceID].positionSets[contigs[sourceID].positionSets.size() - 1][segs[sp].sourceStart + ssp].chromosomeID = targetID;
//			if(fr)
				contigs[sourceID].positionSets[contigs[sourceID].positionSets.size() - 1][segs[sp].sourceStart + ssp].chromosomeOffset = segs[sp].targetStart + ssp;
//			else
//				contigs[sourceID].positionSets[contigs[sourceID].positionSets.size() - 1][segs[sp].sourceStart + ssp].chromosomeOffset = segs[sp].targetStart + segs[sp].size - 1 - ssp;
		}		
//	if(fr == 0 && contigs[sourceID].adjusted == 0)
//	{
//		contigs[sourceID].adjusted = 1;
//		reverse(contigs[sourceID].nucleotides.begin(), contigs[sourceID].nucleotides.end());
//		complement(contigs, sourceID);
//	}
}

void loadContiAli(ifstream & ca, vector<Seq> & contigs, int chromosomeID)
{
	string buf;
	unsigned int targetID, targetStart, targetEnd, targetGap, sourceID, sourceStart, sourceEnd, sourceGap, sourceSize, fr, targetSize;
	int i, cp, pp, ppp, seqID, realSourceID;
	vector<Segment> segs;

	sourceID = -1;
	if(ca.is_open())
	{
		while(ca.good())
		{
			getline(ca, buf);
			if(buf[0] == 0) 
			{
				if(keepPositions(contigs, sourceID, segs, CONTIG_THRESHOLD) == 0)
				{
					contigs[sourceID].positionSets.erase(contigs[sourceID].positionSets.end() - 1);// have to keep or discard the last alignment
				}
				break;
			}

			realSourceID = parseBLAT(buf, targetID, targetStart, targetEnd, targetGap, sourceID, sourceStart, sourceEnd, sourceGap, sourceSize, segs, fr, targetSize);
//			cout << targetID << ", " << targetStart << ", " << targetEnd << ", " << targetGap << ", " << sourceID << ", " << sourceStart << ", " << sourceEnd << ", " << sourceGap << ", " << sourceSize << ", " << segs[0].sourceStart << ", " << segs[0].targetStart << ", " << segs[0].size << ", " << fr << endl;
//			if(segs.size() > 1)
//				cout << segs[1].sourceStart << ", " << segs[1].targetStart << ", " << segs[1].size << endl;

			if((double) (sourceEnd - sourceStart - sourceGap) / sourceSize >= INIT_CONTIG_THRESHOLD && (double) (targetEnd - targetStart - targetGap) / (targetEnd - targetStart) >= INIT_CONTIG_THRESHOLD && sourceSize >= 200)// make the change
//			if((double) (sourceEnd - sourceStart) / sourceSize >= INIT_CONTIG_THRESHOLD)
			{
				updateContig(contigs, sourceID, targetID, segs, fr, CONTIG_THRESHOLD);
			}
		}
	}
	else
	{
		cout << "CANNOT OPEN FILE!" << endl;
		exit(-1);
	}

//	for(cp = 0; cp < contigs.size(); cp ++)
//	{
//		cout << "contig " << contigs[cp].ID << "; ";
//		for(pp = 0; pp < contigs[cp].positionSets.size(); pp ++)
//		{
//			cout << "pos " << contigs[cp].positionSets[pp][0].chromosomeOffset << ", ";
//		}
//		cout << endl;
//	}

//        for(cp = 0; cp < contigs.size(); cp ++)
//                for(pp = 0; pp < contigs[cp].positionSets.size(); pp ++)
//			for(ppp = 0; ppp < contigs[cp].positionSets[pp].size(); ppp ++)
//                	{
//               	        	cout << "<" << contigs[cp].positionSets[pp][ppp].chromosomeID << "," << contigs[cp].positionSets[pp][ppp].chromosomeOffset << ">";
//                        	if((ppp + 1) % 60 == 0 || ppp == contigs[cp].positionSets[pp].size() - 1)
//                                	cout << endl;
//                	}

// output initial contigs
//	ofstream out;
//      string s = "initial_contigs." + itoa(chromosomeID) + ".checked.fa";
//      out.open(s.c_str());
//	seqID = 0;

//      for(cp = 0; cp < contigs.size(); cp ++)
//      {
//		if(contigs[cp].positionSets.size() > 0)
//		{
//              	out << ">" << seqID ++ << endl;
//                	for(pp = 0; pp < contigs[cp].nucleotides.size(); pp ++)
//                	{
//                	        out << contigs[cp].nucleotides[pp];
//                	        if((pp + 1) % 60 == 0 || pp == contigs[cp].nucleotides.size() - 1)
//        	                        out << endl;
//	                }
//		}
//        }
}

void reverseComplement(vector<char> & seq)
{
	int np;

	reverse(seq.begin(), seq.end());
	for(np = 0; np < seq.size(); np ++)
		if(seq[np] == 'A') seq[np] = 'T';
		else if(seq[np] == 'C') seq[np] = 'G';
		else if(seq[np] == 'G') seq[np] = 'C';
		else if(seq[np] == 'T') seq[np] = 'A';
		//do nothing for 'N'
}

void printInitialContigs(ofstream & out, Seq & seq, int seqID)
{
	int np;

	if(seq.outputted == 0)
	{
		out << ">" << seqID << endl;
		for(np = 0; np < seq.nucleotides.size(); np ++)
		{
			out << seq.nucleotides[np];
			if((np + 1) % 60 ==0 || np == seq.nucleotides.size() - 1)
				out << endl;
		}
		seq.outputted = 1;
	}
}

void updateGenomeWithContig(vector<Seq> & seqs, vector<vector<Base> > & genome, int chrID)
{
	int sp, pp, ppp, npp, gp, tag, seqID, size, np;
	unsigned int chromosomeID, chromosomeOffset, nextID, nextOffset, nextItem, cpp, i, previousIDBak, previousOffsetBak, previousItemBak;
	Base b;
	ContiMer p, q;
	KMer k;
	char nucleotide;

	previousIDBak = previousOffsetBak = previousItemBak = -1;
	for(sp = 0; sp < seqs.size(); sp ++)
	{
		pp = 0;
cont:
		size = seqs[sp].positionSets.size() == 0 ? 0 : 1;
		for(; pp < seqs[sp].positionSets.size() ; pp ++)// make the change
		{
			//for(spp = 0; spp < sp; spp ++) // it is too time-consuming to iterate different seqs with different positions for one seq, so function compatible enables merging kmers from with different contig ids
			for(ppp = 0; ppp < pp; ppp ++)
				if(abs((int)(seqs[sp].positionSets[pp][0].chromosomeOffset - seqs[sp].positionSets[ppp][0].chromosomeOffset)) < seqs[sp].nucleotides.size())
				{
					pp ++;
					goto cont;
				}
			for(ppp = 0; ppp < seqs[sp].positionSets[pp].size() - 1; ppp ++)
			{
				if(seqs[sp].positionSets[pp][ppp].chromosomeID != -1)
				{
					chromosomeID = seqs[sp].positionSets[pp][ppp].chromosomeID;
					chromosomeOffset = seqs[sp].positionSets[pp][ppp].chromosomeOffset;
					if(genome[chromosomeID][chromosomeOffset].contiMer.size() >= 2)
					{
						pp ++;
						goto cont;
					}
				}
			}
//			cout << "--------" << pp << "--------" << endl;
			tag = 0;
			if(seqs[sp].frSets[pp] == 1)
			{
				reverseComplement(seqs[sp].nucleotides);
				tag = 1;
			}
//has to always adjust nucleotides according to fr
			seqs[sp].outputted = 1;
//			printInitialContigs(out, seqs[sp], sp);
			for(ppp = 0; ppp < seqs[sp].positionSets[pp].size() - 1; ppp ++)
			{
//				if(seqs[sp].positionSets[pp][ppp].chromosomeOffset == 207630095)
//					cout << "CONTIG TRACED: " << sp << ", " << seqs[sp].positionSets[pp][0].chromosomeOffset  << endl;

				if(seqs[sp].positionSets[pp][ppp].chromosomeID != -1)
				{
					chromosomeID = seqs[sp].positionSets[pp][ppp].chromosomeID;
					chromosomeOffset = seqs[sp].positionSets[pp][ppp].chromosomeOffset;
					nextID = seqs[sp].positionSets[pp][ppp + 1].chromosomeID;
					nextOffset = seqs[sp].positionSets[pp][ppp + 1].chromosomeOffset;
//					nextItem = genome[nextID][nextOffset].nextPosition.size();
					nucleotide = seqs[sp].nucleotides[ppp];
					if(nextID == -1)// insertion to genome
					{
						for(npp = ppp + 2; npp < seqs[sp].positionSets[pp].size(); npp ++)
						{
							if(seqs[sp].positionSets[pp][npp].chromosomeID != -1)
							{
								nextID = seqs[sp].positionSets[pp][npp].chromosomeID;
								nextOffset = seqs[sp].positionSets[pp][npp].chromosomeOffset;
								if(seqs[sp].positionSets[pp][npp].chromosomeID == seqs[sp].positionSets[pp][ppp].chromosomeID && npp - ppp < SI)// small insertion
								{
//									cout << "SI" << endl;
//									for(cpp = chromosomeOffset; cpp < nextOffset; cpp ++)
									{
										//p.previousID = previousIDBak;
										//p.previousOffset = previousOffsetBak;
										//p.previousItem = previousItemBak;
										p.nextID = nextID;
										p.nextOffset = nextOffset;
										p.nextItem = genome[p.nextID][p.nextOffset].contiMer.size();
										p.contigID = sp;
										p.contigOffset = ppp;
										p.nucleotide = nucleotide;
										genome[chromosomeID][chromosomeOffset].contiMer.push_back(p);

										previousIDBak = chromosomeID;
										previousOffsetBak = chromosomeOffset;
										previousItemBak = genome[chromosomeID][chromosomeOffset].contiMer.size() - 1;
										//k.contigOffset = k.contigOffset0 = k.chromosomeID0 = k.chromosomeOffset0 = -1;
										//genome[chromosomeID][chromosomeOffset].kMer.push_back(k);
										//genome[chromosomeID][cpp].nextID = chromosomeID;
										//genome[chromosomeID][cpp].nextOffset = cpp + 1;
									}
//									cout << "IS" << endl;
								}
								else// large insertion
								{
//									cout << "LI" << endl;
									//p.previousID = previousIDBak;
									//p.previousOffset = previousOffsetBak;
									//p.previousItem = previousItemBak;
									p.nextID = chromosomeID;
									p.nextOffset = genome[chromosomeID].size();
									p.nextItem = 0;
									p.contigID = sp;
									p.contigOffset = ppp;
									p.nucleotide = nucleotide;
									genome[chromosomeID][chromosomeOffset].contiMer.push_back(p);
									previousIDBak = chromosomeID;
									previousOffsetBak = chromosomeOffset;
									previousItemBak = genome[chromosomeID][chromosomeOffset].contiMer.size() - 1;
									//k.contigOffset = k.contigOffset0 = k.chromosomeID0 = k.chromosomeOffset0 = -1;
									//genome[chromosomeID][chromosomeOffset].kMer.push_back(k);
									//genome[chromosomeID][chromosomeOffset].nextID = chromosomeID;
									//genome[chromosomeID][chromosomeOffset].nextOffset = genome[chromosomeID].size();
									for(i = 0; i < npp - ppp - 2; i ++)
									{
										b.contiMer.clear();
										b.nucleotide = seqs[sp].nucleotides[ppp + 1 + i];
										//p.previousID = previousIDBak;
										//p.previousOffset = previousOffsetBak;
										//p.previousItem = previousItemBak;
										p.nextID = chromosomeID;
										p.nextOffset = genome[chromosomeID].size() + 1;
										p.nextItem = 0;
										p.contigID = sp;
										p.contigOffset = ppp + 1 + i;
										p.nucleotide = b.nucleotide;
										b.contiMer.push_back(p);

										previousIDBak = chromosomeID;
										previousOffsetBak = genome[chromosomeID].size();
										previousItemBak = 0;
										//k.contigOffset = k.contigOffset0 = k.chromosomeID0 = k.chromosomeOffset0 = -1;
										//b.kMer.push_back(k);
										//b.nextID = chromosomeID;
										//b.nextOffset = genome[chromosomeID].size();
										genome[chromosomeID].push_back(b);
									}
									b.contiMer.clear();
									b.nucleotide = seqs[sp].nucleotides[npp - 1];
									//p.previousID = previousIDBak;
									//p.previousOffset = previousOffsetBak;
									//p.previousItem = previousItemBak;
									p.nextID = nextID; //contigs[cp].positionSets[pp][npp].chromosomeID;
									p.nextOffset = nextOffset; //contigs[cp].positionSets[pp][npp].chromosomeOffset;
									p.nextItem = genome[p.nextID][p.nextOffset].contiMer.size();
									p.contigID = sp;
									p.contigOffset = npp - 1;
									p.nucleotide = b.nucleotide;
									b.contiMer.push_back(p);

									previousIDBak = chromosomeID;
									previousOffsetBak = genome[chromosomeID].size();
									previousItemBak = 0;
									//k.contigOffset = k.contigOffset0 = k.chromosomeID0 = k.chromosomeOffset0 = -1;
									//b.kMer.push_back(k);
									genome[chromosomeID].push_back(b);
									//genome[chromosomeID][genome[chromosomeID].size() - 1].nextID = nextID;
									//genome[chromosomeID][genome[chromosomeID].size() - 1].nextOffset = nextOffset;
//									cout << "IL" << endl;
								}
								ppp = npp - 1;
								break;
							}
						}
					}
					else if(nextID == chromosomeID && nextOffset != chromosomeOffset + 1)// deletion from genome
					{
						if(nextOffset - chromosomeOffset < SD)// small deletion
						{
//							cout << "SD" << endl;
							for(cpp = chromosomeOffset; cpp < nextOffset; cpp ++)
							{
								//p.previousID = previousIDBak;
								//p.previousOffset = previousOffsetBak;
								//p.previousItem = previousItemBak;
								p.nextID = chromosomeID;
								p.nextOffset = cpp + 1;
								p.nextItem = genome[p.nextID][p.nextOffset].contiMer.size();
								p.contigID = sp;
								p.contigOffset = ppp;
								if(cpp == chromosomeOffset) p.nucleotide = nucleotide;
								else p.nucleotide = genome[chromosomeID][cpp].nucleotide;
								genome[chromosomeID][cpp].contiMer.push_back(p);

								previousIDBak = chromosomeID;
								previousOffsetBak = cpp;
								previousItemBak = genome[chromosomeID][cpp].contiMer.size() - 1;
								//k.contigOffset = k.contigOffset0 = k.chromosomeID0 = k.chromosomeOffset0 = -1;
								//genome[chromosomeID][cpp].kMer.push_back(k);
								//genome[chromosomeID][cpp].nextID = chromosomeID;
								//genome[chromosomeID][cpp].nextOffset = cpp + 1;
							}
//							cout << "DS" << endl;
						}
						else// large deletion
						{
//							cout << "LD" << endl;
							//p.previousID = previousIDBak;
							//p.previousOffset = previousOffsetBak;
							//p.previousItem = previousItemBak;
							p.nextID = nextID;
							p.nextOffset = nextOffset;
							p.nextItem = genome[p.nextID][p.nextOffset].contiMer.size();
							p.contigID = sp;
							p.contigOffset = ppp;
							p.nucleotide = nucleotide;
							genome[chromosomeID][chromosomeOffset].contiMer.push_back(p);

							previousIDBak = chromosomeID;
							previousOffsetBak = chromosomeOffset;
							previousItemBak = genome[chromosomeID][chromosomeOffset].contiMer.size() - 1;
							//k.contigOffset = k.contigOffset0 = k.chromosomeID0 = k.chromosomeOffset0 = -1;
							//genome[chromosomeID][chromosomeOffset].kMer.push_back(k);
							//genome[chromosomeID][chromosomeOffset].nextID = nextID;
							//genome[chromosomeID][chromosomeOffset].nextOffset = nextOffset;
//							cout << "DL" << endl;
						}
					}
					else// ordinary case
					{
						//p.previousID = previousIDBak;
						//p.previousOffset = previousOffsetBak;
						//p.previousItem = previousItemBak;
						p.nextID = nextID;
						p.nextOffset = nextOffset;
						p.nextItem = genome[p.nextID][p.nextOffset].contiMer.size();
						p.contigID = sp;
						p.contigOffset = ppp;
						p.nucleotide = nucleotide;
						genome[chromosomeID][chromosomeOffset].contiMer.push_back(p);
						previousIDBak = chromosomeID;
						previousOffsetBak = chromosomeOffset;
						previousItemBak = genome[chromosomeID][chromosomeOffset].contiMer.size() - 1;
						//k.contigOffset = k.contigOffset0 = k.chromosomeID0 = k.chromosomeOffset0 = -1;
						//genome[chromosomeID][chromosomeOffset].kMer.push_back(k);
						//genome[chromosomeID][chromosomeOffset].nextID = nextID;
						//genome[chromosomeID][chromosomeOffset].nextOffset = nextOffset;
					}
				}
			}
			if(nextID != -1)
			{
				//p.previousID = previousIDBak;
				//p.previousOffset = previousOffsetBak;
				//p.previousItem = previousItemBak;
				p.nextID = p.nextOffset = p.nextItem = -1;
				p.contigID = sp;
				p.contigOffset = ppp;
				p.nucleotide = genome[nextID][nextOffset].nucleotide;
				genome[nextID][nextOffset].contiMer.push_back(p);

				//k.contigOffset = k.contigOffset0 = k.chromosomeID0 = k.chromosomeOffset0 = -1;
				//genome[nextID][nextOffset].kMer.push_back(k);
			}
			else// chromosomeID != -1 and nextID == -1
			{
				//p.previousID = previousIDBak;
				//p.previousOffset = previousOffsetBak;
				//p.previousItem = previousItemBak;
				p.nextID = p.nextOffset = p.nextItem = -1;
				p.contigID = sp;
				p.contigOffset = ppp;
				p.nucleotide = genome[chromosomeID][chromosomeOffset].nucleotide;
				genome[chromosomeID][chromosomeOffset].contiMer.push_back(p);

				//k.contigOffset = k.contigOffset0 = k.chromosomeID0 = k.chromosomeOffset0 = -1;
				//genome[chromosomeID][chromosomeOffset].kMer.push_back(k);
			}
			if(tag == 1) reverseComplement(seqs[sp].nucleotides);
		}

//		cout << "*************************" << endl;
//		for(gp = 0; gp < genome.size(); gp ++)
//		{
//			for(cpp = 0; cpp < genome[gp].size(); cpp ++)
//			{
//				cout << "[";
//				for(pp = 0; pp < genome[gp][cpp].contiMer.size(); pp ++)
//					cout << "<" << genome[gp][cpp].contiMer[pp].contigOffset << "| " << genome[gp][cpp].contiMer[pp].previousID << ", " << genome[gp][cpp].contiMer[pp].previousOffset << ", " << genome[gp][cpp].contiMer[pp].previousItem << "| " << genome[gp][cpp].contiMer[pp].nextID << ", " << genome[gp][cpp].contiMer[pp].nextOffset << ", " << genome[gp][cpp].contiMer[pp].nextItem << ">";
//				cout << "]";
//			}
//			cout << endl;
//		}
//        	cout << "-------------------------" << endl;
//        	for(gp = 0; gp < genome.size(); gp ++)
//        	{
//                	for(cpp = 0; cpp < genome[gp].size(); cpp ++)
//                	{
//                        	cout << "[";
//                        	for(pp = 0; pp < genome[gp][cpp].kMer.size(); pp ++)
//                                	cout << "<" << genome[gp][cpp].kMer[pp].contigOffset << ", " << genome[gp][cpp].kMer[pp].contigOffset0 << ", " << genome[gp][cpp].kMer[pp].chromosomeID0 << ", " << genome[gp][cpp].kMer[pp].chromosomeOffset0 <<  "| " << genome[gp][cpp].kMer[pp].previous[0].previousID << ", " << genome[gp][cpp].kMer[pp].previous[0].previousOffset << ", " << genome[gp][cpp].kMer[pp].previous[0].previousItem << "| " << genome[gp][cpp].kMer[pp].next[0].nextID << ", " << genome[gp][cpp].kMer[pp].next[0].nextOffset << ", " << genome[gp][cpp].kMer[pp].next[0].nextItem << ">";
//                        	cout << "]";
//                	}
//        	}
//		cout << endl;
//		cout << "*************************" << endl;
	}

        ofstream out;
        string s = "tmp/_initial_contigs." + itoa(chrID) + ".fa";
        out.open(s.c_str());
	vector<vector<char> > contigs;
	vector<char> contig;
	vector<int> outputted;
	vector<int> maxOutputted;
	int cp, frag;

	int IDBak = -1;
        for(sp = 0; sp < seqs.size(); sp ++)
        {
		if(seqs[sp].ID != IDBak)
		{
			maxOutputted.push_back(0);
			outputted.push_back(0);
			contigs.push_back(contig);
			IDBak = seqs[sp].ID;
		}

		maxOutputted[maxOutputted.size() - 1] ++;
		outputted[outputted.size() - 1] = outputted[outputted.size() - 1] + seqs[sp].outputted;
		for(np = 0; np < seqs[sp].nucleotides.size(); np ++)
			contigs[contigs.size() - 1].push_back(seqs[sp].nucleotides[np]);
	}

	for(cp = 0; cp < contigs.size(); cp ++)
	{
		if((double) outputted[cp] / (double) maxOutputted[cp] >= CONTIG_THRESHOLD)
		{
	                out << ">" << cp << endl;
	                for(np = 0; np < contigs[cp].size(); np ++)
	                {
	                        out << contigs[cp][np];
	                        if((np + 1) % 60 == 0 || np == contigs[cp].size() - 1)
	                                out << endl;
	                }
		}
        }
}

void loadContigAlignment(vector<vector<Base> > & genome, int chromosomeID)
{
        vector<Seq> contigs;
	ifstream c, ca;
	string s;

	c.open("tmp/_contigs.fa");
	s = "tmp/_contigs_genome." + itoa(chromosomeID) + ".psl";
	ca.open(s.c_str());
        loadSeq(c, contigs);
        loadContiAli(ca, contigs, chromosomeID);
        updateGenomeWithContig(contigs, genome, chromosomeID);
}

int loadReadAli(ifstream & ra, vector<Seq> & reads, int & aliStartID, int & seqStartID)
{
	string buf;
	unsigned int targetID1, targetStart1, targetEnd1, targetGap1, sourceID1, sourceStart1, sourceEnd1, sourceGap1, sourceSize1, fr1, targetID2, targetStart2, targetEnd2, targetGap2, sourceID2, sourceStart2, sourceEnd2, sourceGap2, sourceSize2, fr2;
	vector<Segment> segs1, segs2;
	int gp, cpp, pp;
	int count = 0;

	sourceIDBak = -1;
	if(ra.is_open())
	{
		while(ra.good())
		{
			getline(ra, buf);
			if(buf[0] == 0) break;
			if(buf[0] == '@') continue;
			//segs1.clear();// it is important to push_back first and then clear
			parseBOWTIE(buf, targetID1, targetStart1, targetEnd1, targetGap1, sourceID1, sourceStart1, sourceEnd1, sourceGap1, sourceSize1, segs1, fr1);
//			cout << targetID1 << ", " << targetStart1 << ", " << targetEnd1 << ", " << targetGap1 << ", " << sourceID1 << ", " << sourceStart1 << ", " << sourceEnd1 << ", " << sourceGap1 << ", " << sourceSize1 << ", " << segs1[0].sourceStart << ", " << segs1[0].targetStart << ", " << segs1[0].size << ", " << fr1 << endl;
			getline(ra, buf);
			if(buf[0] == 0)
			{
				cout << "BROKEN BOWTIE FILE" << endl;
				exit(-1);
			}
			//segs2.clear();
			parseBOWTIE(buf, targetID2, targetStart2, targetEnd2, targetGap2, sourceID2, sourceStart2, sourceEnd2, sourceGap2, sourceSize2, segs2, fr2);
//			cout << targetID2 << ", " << targetStart2 << ", " << targetEnd2 << ", " << targetGap2 << ", " << sourceID2 << ", " << sourceStart2 << ", " << sourceEnd2 << ", " << sourceGap2 << ", " << sourceSize2 << ", " << segs2[0].sourceStart << ", " << segs2[0].targetStart << ", " << segs2[0].size << ", " << fr2 << endl;

			if(sourceID1 < aliStartID) {/*cout << "continued" << endl;*/ continue;}
			if(sourceID1 > seqStartID) {/*cout << "unfinished" << endl;*/ return 0;}

			if(targetID1 != -1 && targetID2 != -1 && (double) (sourceEnd1 - sourceStart1 - sourceGap1) / /*(sourceEnd1 - sourceStart1)*/ sourceSize1 >= THRESHOLD && (double) (targetEnd1 - targetStart1 - targetGap1) / (targetEnd1 - targetStart1) >= THRESHOLD && (double) (sourceEnd2 - sourceStart2 - sourceGap2) / /*(sourceEnd2 - sourceStart2)*/ sourceSize2 >= THRESHOLD && (double) (targetEnd2 - targetStart2 - targetGap2) / (targetEnd2 - targetStart2) >= THRESHOLD)// 0.8 to 0 for test
			{
				updateContig(reads, (sourceID1 - aliStartID) * 2, targetID1, segs1, fr1, THRESHOLD);
				updateContig(reads, (sourceID2 - aliStartID) * 2 + 1, targetID2, segs2, fr2, THRESHOLD);
				count ++;
			}
			segs1.clear();
			segs2.clear();
		}
		return 1;
	}
	else
	{
		cout << "CANNOT OPEN FILE!" << endl;
		exit(-1);
	}
}

int getDiff(vector<char> s1, vector<char> s2)
{
	int diff = 0, i;

	if(s1.size() == 0 || s2.size() == 0) return -1;

	for(i = 0; i < s1.size(); i ++)
		if(s1[i] != s2[i])
			diff ++;
	return diff;
}

//allowed pair distance < distanceHigh - distanceLow so that PE de Bruijn graph's branch solving function works
//distanceHigh - distanceLow > 50
int compatible(KMer & k1, KMer & k2, int insertVariation, int mrl)
{
//	cout << endl << "---------KMER COMPARE---------" << endl;
//	cout << "contigID: " << k1.contigID << " vs " << k2.contigID << endl;
//	cout << "contigOffset: " << k1.contigOffset << " vs " << k2.contigOffset << ": " << abs((int)(k1.contigOffset - k2.contigOffset)) << endl;
//	cout << "contigID0: " << k1.contigID0 << " vs " << k2.contigID0 << endl;
//	cout << "contigOffset0 " << k1.contigOffset0 << " vs " << k2.contigOffset0 << ": " << abs((int)(k1.contigOffset0 - k2.contigOffset0)) << endl;
//	cout << "chromosomeID0: " << k1.chromosomeID0 << " vs " << k2.chromosomeID0 << endl;
//	cout << "chromosomeOffset0: " << k1.chromosomeOffset0 << " vs " << k2.chromosomeOffset0 << ": " << abs((int)(k1.chromosomeOffset0 - k2.chromosomeOffset0)) << endl;

	if(
	(k1.contigID == -1 || k2.contigID == -1 || (k1.contigID != -1 && k1.contigID == k2.contigID && abs((int)(k1.contigOffset - k2.contigOffset)) <= 5 * EP) 
#ifdef OPTIMIZATION
	|| (k1.contigID != -1 && k2.contigID != -1 && k1.contigID != k2.contigID)
#endif
	) &&
// it is hard to judge if two overlapping contigs should be joined or are just repetitive. Here I find it performs better to join them.
	(k1.contigID0 == -1 || k2.contigID0 == -1 || (k1.contigID0 != -1 && k1.contigID0 == k2.contigID0 && abs((int)(k1.contigOffset0 - k2.contigOffset0)) <= 2 * insertVariation + 5 * EP) 
#ifdef OPTIMIZATION
	|| (k1.contigID0 != -1 && k2.contigID0 != -1 && k1.contigID0 != k2.contigID0)
#endif
	) &&
	(k1.chromosomeID0 == -1 || k2.chromosomeID0 == -1 || (k1.chromosomeID0 != -1 && k1.chromosomeID0 == k2.chromosomeID0 && abs((int)(k1.chromosomeOffset0 - k2.chromosomeOffset0)) <= 2 * insertVariation + 5 * EP))
	)
	{
//		if(k2.nextID == -1)
//		{
//			k2.nextID = k1.nextID;
//			k2.nextOffset = k1.nextOffset;
//			k2.nextItem = k1.nextItem;
//		}
//		if(k2.previousID == -1)
//		{
//			k2.previousID = k1.previousID;
//			k2.previousOffset = k1.previousOffset;
//			k2.previousItem = k1.previousItem;
//		}
//		cout << "compatible" << endl;
//		cout << "------------------------" << endl;
		return 1;
	}
	else
	{
//		cout << "not compatible" << endl;
//		cout << "------------------------" << endl;
		return 0;
	}
}

int nextCompatible(Next next1, Next next2)
{
	if(next1.nextID == -1 || next2.nextID == -1)
	{
		cout << "KMER ERROR" << endl;
		exit(-1);
	}

	if(next1.nextID == next2.nextID && next1.nextOffset == next2.nextOffset && next1.nextItem == next2.nextItem)
		return 1;
	return 0;
}

int previousCompatible(Previous previous1, Previous previous2)
{
        if(previous1.previousID == -1 || previous2.previousID == -1)
        {
                cout << "KMER ERROR" << endl;
                exit(-1);
        }

        if(previous1.previousID == previous2.previousID && previous1.previousOffset == previous2.previousOffset && previous1.previousItem == previous2.previousItem)
                return 1;
        return 0;
}

void updateKBases(KMer k, KMer & k0)
{
//	int np, size;

//	size = k.s.size() < 1 ? k.s.size() : 1;
//	for(np = 0; np < size/*k.s.size()*/; np ++)
//	{
//		if(np >= k0.A.size())
//		{
//			k0.A.push_back(0x00);
//			k0.C.push_back(0x00);
//			k0.G.push_back(0x00);
//			k0.T.push_back(0x00);
//			k0.N.push_back(0x00);
//		}

//		switch(k.s[np])
//		{
//			case 'A': if(k0.A[np] != 0xff) k0.A[np] ++; break;
//			case 'C': if(k0.C[np] != 0xff) k0.C[np] ++; break;
//			case 'G': if(k0.G[np] != 0xff) k0.G[np] ++; break;
//			case 'T': if(k0.T[np] != 0xff) k0.T[np] ++; break;
//			default: if(k0.N[np] != 0xff) k0.N[np] ++;
//		}
//	}

	if(k.s.size() > 0)
	        switch(k.s[0])
	        {
	                case 'A': k0.A ++; break;
	                case 'C': k0.C ++; break;
	                case 'G': k0.G ++; break;
	                case 'T': k0.T ++; break;
	                default: k0.N ++;
	        }
}

void updateKMer(vector<vector<Base> > & genome, unsigned int chromosomeID, unsigned int chromosomeOffset, unsigned int nextID, unsigned int nextOffset, unsigned int chromosomeID0, unsigned int chromosomeOffset0, unsigned int nextID0, unsigned int nextOffset0, vector<char> s, vector<char> nextS, int insertVariation, int mrl)
{
	int ip, ipp, ip0, np, pp, cip, nip, size;
	KMer k1, k2;
	ContiMer p;
	vector<unsigned int> nextItem, chromosomeItem;
	Next next;
	Previous previous;

	k1.traversed = 0;
	k1.s = s;
	k1.chromosomeID0 = chromosomeID0;
	k1.chromosomeOffset0 = chromosomeOffset0;
//	size = s.size() < 1 ? s.size() : 1;
//	for(np = 0; np < size/*s.size()*/; np ++)
//	{
//		k1.A.push_back(0x00);
//		k1.C.push_back(0x00);
//		k1.G.push_back(0x00);
//		k1.T.push_back(0x00);
//		k1.N.push_back(0x00);
//	}
	k1.A = k1.C = k1.G = k1.T = k1.N = 0;
//	k1.nextID = -1;//nextID;
//	k1.nextOffset = -1;//nextOffset;
//	k1.nextItem = -1;//k1.nextID == -1 ? -1 : genome[k1.nextID][k1.nextOffset].kMer.size();
//	k1.previousID = k1.previousOffset = k1.previousItem = -1;
	k1.coverage = 1;

	if(genome[chromosomeID][chromosomeOffset].contiMer.size() == 0 && (chromosomeID0 == -1 || genome[chromosomeID0][chromosomeOffset0].contiMer.size() == 0))
	{
		k1.contigID = -1;
		k1.contigID0 = -1;
		k1.contigOffset = -1;
		k1.contigOffset0 = -1;
		for(ipp = 0; ipp < genome[chromosomeID][chromosomeOffset].kMer.size(); ipp ++)
			if(compatible(k1, genome[chromosomeID][chromosomeOffset].kMer[ipp], insertVariation, mrl))
				break;
		if(ipp == genome[chromosomeID][chromosomeOffset].kMer.size())
		{
//			cout << "insert to " << chromosomeID << ", " << chromosomeOffset << ", " << ipp << endl;
			genome[chromosomeID][chromosomeOffset].kMer.push_back(k1);
			updateKBases(k1, genome[chromosomeID][chromosomeOffset].kMer[ipp]);
		}
		else
		{
			genome[chromosomeID][chromosomeOffset].kMer[ipp].coverage ++;
//			cout << "update " << chromosomeID << ", " << chromosomeOffset << ", " << ipp << endl;
			updateKBases(k1, genome[chromosomeID][chromosomeOffset].kMer[ipp]);
		}
//		chromosomeItem = ipp;
		chromosomeItem.push_back(ipp);
		goto cont;
	}

	if(genome[chromosomeID][chromosomeOffset].contiMer.size() != 0 && (chromosomeID0 == -1 || genome[chromosomeID0][chromosomeOffset0].contiMer.size() == 0))
	{
		k1.contigID0 = -1;
		k1.contigOffset0 = -1;
		for(ip = 0; ip < genome[chromosomeID][chromosomeOffset].contiMer.size(); ip ++)
		{
			k1.contigID = genome[chromosomeID][chromosomeOffset].contiMer[ip].contigID;
			k1.contigOffset = genome[chromosomeID][chromosomeOffset].contiMer[ip].contigOffset;
			for(ipp = 0; ipp < genome[chromosomeID][chromosomeOffset].kMer.size(); ipp ++)
				if(compatible(k1, genome[chromosomeID][chromosomeOffset].kMer[ipp], insertVariation, mrl))
					break;
			if(ipp == genome[chromosomeID][chromosomeOffset].kMer.size())
			{
//				cout << "insert to " << chromosomeID << ", " << chromosomeOffset << ", " << ipp << endl;
				genome[chromosomeID][chromosomeOffset].kMer.push_back(k1);
				updateKBases(k1, genome[chromosomeID][chromosomeOffset].kMer[ipp]);
			}
			else
			{
				genome[chromosomeID][chromosomeOffset].kMer[ipp].coverage ++;
//				cout << "update " << chromosomeID << ", " << chromosomeOffset << ", " << ipp << endl;
				updateKBases(k1, genome[chromosomeID][chromosomeOffset].kMer[ipp]);
			}
			chromosomeItem.push_back(ipp);
		}
//		chromosomeItem = ipp;
		goto cont;
	}

	if(genome[chromosomeID][chromosomeOffset].contiMer.size() == 0 && (chromosomeID0 != -1 && genome[chromosomeID0][chromosomeOffset0].contiMer.size() != 0))
	{
		k1.contigID = -1;
		k1.contigOffset = -1;
		for(ip0 = 0; ip0 < genome[chromosomeID0][chromosomeOffset0].contiMer.size(); ip0 ++)
		{
			k1.contigID0 = genome[chromosomeID0][chromosomeOffset0].contiMer[ip0].contigID;
			k1.contigOffset0 = genome[chromosomeID0][chromosomeOffset0].contiMer[ip0].contigOffset;
			for(ipp = 0; ipp < genome[chromosomeID][chromosomeOffset].kMer.size(); ipp ++)
				if(compatible(k1, genome[chromosomeID][chromosomeOffset].kMer[ipp], insertVariation, mrl))
					break;
			if(ipp == genome[chromosomeID][chromosomeOffset].kMer.size())
			{
//				cout << "insert to " << chromosomeID << ", " << chromosomeOffset << ", " << ipp << endl;
				genome[chromosomeID][chromosomeOffset].kMer.push_back(k1);
				updateKBases(k1, genome[chromosomeID][chromosomeOffset].kMer[ipp]);
			}
			else
			{
				genome[chromosomeID][chromosomeOffset].kMer[ipp].coverage ++;
//				cout << "update " << chromosomeID << ", " << chromosomeOffset << ", " << ipp << endl;
				updateKBases(k1, genome[chromosomeID][chromosomeOffset].kMer[ipp]);
			}
			chromosomeItem.push_back(ipp);
		}
//		chromosomeItem = ipp;
		goto cont;
	}

	if(genome[chromosomeID][chromosomeOffset].contiMer.size() != 0 && (chromosomeID0 != -1 && genome[chromosomeID0][chromosomeOffset0].contiMer.size() != 0))
	{
		for(ip = 0; ip < genome[chromosomeID][chromosomeOffset].contiMer.size(); ip ++)
			for(ip0 = 0; ip0 < genome[chromosomeID0][chromosomeOffset0].contiMer.size(); ip0 ++)
			{
				k1.contigID = genome[chromosomeID][chromosomeOffset].contiMer[ip].contigID;
				k1.contigID0 = genome[chromosomeID0][chromosomeOffset0].contiMer[ip0].contigID;
				k1.contigOffset = genome[chromosomeID][chromosomeOffset].contiMer[ip].contigOffset;
				k1.contigOffset0 = genome[chromosomeID0][chromosomeOffset0].contiMer[ip0].contigOffset;
				for(ipp = 0; ipp < genome[chromosomeID][chromosomeOffset].kMer.size(); ipp ++)
					if(compatible(k1, genome[chromosomeID][chromosomeOffset].kMer[ipp], insertVariation, mrl))
						break;
				if(ipp == genome[chromosomeID][chromosomeOffset].kMer.size())
				{
//					cout << "insert to " << chromosomeID << ", " << chromosomeOffset << ", " << ipp << endl;
					genome[chromosomeID][chromosomeOffset].kMer.push_back(k1);
					updateKBases(k1, genome[chromosomeID][chromosomeOffset].kMer[ipp]);
				}
				else
				{
					genome[chromosomeID][chromosomeOffset].kMer[ipp].coverage ++;
//					cout << "update " << chromosomeID << ", " << chromosomeOffset << ", " << ipp << endl;
					updateKBases(k1, genome[chromosomeID][chromosomeOffset].kMer[ipp]);
				}
				chromosomeItem.push_back(ipp);
			}
//		chromosomeItem = ipp;
		goto cont;
	}

cont:
	k2.traversed = 0;
	k2.s = nextS;
        k2.chromosomeID0 = nextID0;
        k2.chromosomeOffset0 = nextOffset0;
//      k2.previousID = -1;//chromosomeID;
//      k2.previousOffset = -1;//chromosomeOffset;
//      k2.previousItem = -1;//genome[k2.previousID][k2.previousOffset].kMer.size() - 1;
//      k2.nextID = k2.nextOffset = k2.nextItem = -1;
	k2.coverage = 0;
//	size = nextS.size() < 1 ? nextS.size() : 1;
//        for(np = 0; np < size/*nextS.size()*/; np ++)
//        {
//                k2.A.push_back(0x00);
//                k2.C.push_back(0x00);
//                k2.G.push_back(0x00);
//                k2.T.push_back(0x00);
//                k2.N.push_back(0x00);
//        }
	k2.A = k2.C = k2.G = k2.T = k2.N = 0;

        if(genome[nextID][nextOffset].contiMer.size() == 0 && (nextID0 == -1 || genome[nextID0][nextOffset0].contiMer.size() == 0))
        {
		k2.contigID = -1;
		k2.contigID0 = -1;
                k2.contigOffset = -1;
                k2.contigOffset0 = -1;
                for(ipp = 0; ipp < genome[nextID][nextOffset].kMer.size(); ipp ++)
                        if(compatible(k2, genome[nextID][nextOffset].kMer[ipp], insertVariation, mrl))
                                break;
                if(ipp == genome[nextID][nextOffset].kMer.size())
                {
//                      cout << "insert to " << nextID << ", " << nextOffset << ", " << ipp << endl;
                        genome[nextID][nextOffset].kMer.push_back(k2);
                }
                else
		{
//			genome[chromosomeID][chromosomeOffset].kMer[ipp].coverage ++;
//                      cout << "update " << nextID << ", " << nextOffset << ", " << ipp << endl;
		}
//		nextItem = ipp;
		nextItem.push_back(ipp);
		goto conti;
        }

        if(genome[nextID][nextOffset].contiMer.size() != 0 && (nextID0 == -1 || genome[nextID0][nextOffset0].contiMer.size() == 0))
        {
		k2.contigID0 = -1;
                k2.contigOffset0 = -1;
                for(ip = 0; ip < genome[nextID][nextOffset].contiMer.size(); ip ++)
                {
			k2.contigID = genome[nextID][nextOffset].contiMer[ip].contigID;
                        k2.contigOffset = genome[nextID][nextOffset].contiMer[ip].contigOffset;
                        for(ipp = 0; ipp < genome[nextID][nextOffset].kMer.size(); ipp ++)
                                if(compatible(k2, genome[nextID][nextOffset].kMer[ipp], insertVariation, mrl))
                                        break;
                        if(ipp == genome[nextID][nextOffset].kMer.size())
                        {
//                              cout << "insert to " << nextID << ", " << nextOffset << ", " << ipp << endl;
                                genome[nextID][nextOffset].kMer.push_back(k2);
                        }
                        else
			{
//				genome[chromosomeID][chromosomeOffset].kMer[ipp].coverage ++;
//                              cout << "update " << nextID << ", " << nextOffset << ", " << ipp << endl;
			}
			nextItem.push_back(ipp);
                }
//		nextItem = ipp;
		goto conti;
        }

        if(genome[nextID][nextOffset].contiMer.size() == 0 && (nextID0 != -1 && genome[nextID0][nextOffset0].contiMer.size() != 0))
        {
		k2.contigID = -1;
                k2.contigOffset = -1;
                for(ip0 = 0; ip0 < genome[nextID0][nextOffset0].contiMer.size(); ip0 ++)
                {
			k2.contigID0 = genome[nextID0][nextOffset0].contiMer[ip0].contigID;
                        k2.contigOffset0 = genome[nextID0][nextOffset0].contiMer[ip0].contigOffset;
                        for(ipp = 0; ipp < genome[nextID][nextOffset].kMer.size(); ipp ++)
                                if(compatible(k2, genome[nextID][nextOffset].kMer[ipp], insertVariation, mrl))
                                        break;
                        if(ipp == genome[nextID][nextOffset].kMer.size())
                        {
//                              cout << "insert to " << nextID << ", " << nextOffset << ", " << ipp << endl;
                                genome[nextID][nextOffset].kMer.push_back(k2);
                        }
                        else
			{
//				genome[chromosomeID][chromosomeOffset].kMer[ipp].coverage ++;
//                              cout << "update " << nextID << ", " << nextOffset << ", " << ipp << endl;
			}
			nextItem.push_back(ipp);
                }
//		nextItem = ipp;
		goto conti;
        }

        if(genome[nextID][nextOffset].contiMer.size() != 0 && (nextID0 != -1 && genome[nextID0][nextOffset0].contiMer.size() != 0))
	{
                for(ip = 0; ip < genome[nextID][nextOffset].contiMer.size(); ip ++)
                        for(ip0 = 0; ip0 < genome[nextID0][nextOffset0].contiMer.size(); ip0 ++)
                        {
				k2.contigID = genome[nextID][nextOffset].contiMer[ip].contigID;
				k2.contigID0 = genome[nextID0][nextOffset0].contiMer[ip0].contigID;
                                k2.contigOffset = genome[nextID][nextOffset].contiMer[ip].contigOffset;
                                k2.contigOffset0 = genome[nextID0][nextOffset0].contiMer[ip0].contigOffset;
                                for(ipp = 0; ipp < genome[nextID][nextOffset].kMer.size(); ipp ++)
                                        if(compatible(k2, genome[nextID][nextOffset].kMer[ipp], insertVariation, mrl))
                                                break;
                                if(ipp == genome[nextID][nextOffset].kMer.size())
                                {
//                                    cout << "insert to " << nextID << ", " << nextOffset << ", " << ipp << endl;
                                        genome[nextID][nextOffset].kMer.push_back(k2);
                                }
                                else
				{
//					genome[chromosomeID][chromosomeOffset].kMer[ipp].coverage ++;
//                                    cout << "update " << nextID << ", " << nextOffset << ", " << ipp << endl;
				}
				nextItem.push_back(ipp);
                        }
//		nextItem = ipp;
		goto conti;
	}
//it is important to consider the case that one read geneartes two or more k-mers
conti:
	for(cip = 0; cip < chromosomeItem.size(); cip ++)
		for(nip = 0; nip < nextItem.size(); nip ++)
		{
			next.nextID = nextID;
			next.nextOffset = nextOffset;
			next.nextItem = nextItem[nip];
			for(np = 0; np < genome[chromosomeID][chromosomeOffset].kMer[chromosomeItem[cip]].next.size(); np ++)
				if(nextCompatible(next, genome[chromosomeID][chromosomeOffset].kMer[chromosomeItem[cip]].next[np]))
					break;
//decide connectivity between m new k1's and n new k2's
			if(np == genome[chromosomeID][chromosomeOffset].kMer[chromosomeItem[cip]].next.size() && 

			(genome[nextID][nextOffset].kMer[nextItem[nip]].contigID == -1 || 
			genome[chromosomeID][chromosomeOffset].kMer[chromosomeItem[cip]].contigID == -1 || 
			(genome[nextID][nextOffset].kMer[nextItem[nip]].contigID != -1 && genome[nextID][nextOffset].kMer[nextItem[nip]].contigID == genome[chromosomeID][chromosomeOffset].kMer[chromosomeItem[cip]].contigID && abs((int)(genome[nextID][nextOffset].kMer[nextItem[nip]].contigOffset - genome[chromosomeID][chromosomeOffset].kMer[chromosomeItem[cip]].contigOffset)) <= 5 * EP) 
#ifdef OPTIMIZATION
			|| genome[nextID][nextOffset].kMer[nextItem[nip]].contigID != -1 && genome[chromosomeID][chromosomeOffset].kMer[chromosomeItem[cip]].contigID != -1 && genome[nextID][nextOffset].kMer[nextItem[nip]].contigID != genome[chromosomeID][chromosomeOffset].kMer[chromosomeItem[cip]].contigID 
#endif
			) &&
                        (genome[nextID][nextOffset].kMer[nextItem[nip]].contigID0 == -1 ||
                        genome[chromosomeID][chromosomeOffset].kMer[chromosomeItem[cip]].contigID0 == -1 ||
                        (genome[nextID][nextOffset].kMer[nextItem[nip]].contigID0 != -1 && genome[nextID][nextOffset].kMer[nextItem[nip]].contigID0 == genome[chromosomeID][chromosomeOffset].kMer[chromosomeItem[cip]].contigID0 && abs((int)(genome[nextID][nextOffset].kMer[nextItem[nip]].contigOffset0 - genome[chromosomeID][chromosomeOffset].kMer[chromosomeItem[cip]].contigOffset0)) <= 2 * insertVariation + 5 * EP) 
#ifdef OPTIMIZATION
			|| genome[nextID][nextOffset].kMer[nextItem[nip]].contigID0 != -1 && genome[chromosomeID][chromosomeOffset].kMer[chromosomeItem[cip]].contigID0 != -1 && genome[nextID][nextOffset].kMer[nextItem[nip]].contigID0 != genome[chromosomeID][chromosomeOffset].kMer[chromosomeItem[cip]].contigID0
#endif
			))

			{
//				cout << "connect [" << chromosomeID << ", " << chromosomeOffset << ", " << cip << ", " << np << "] to [" << nextID << ", " << nextOffset << ", " << nip << "]" << endl;
				genome[chromosomeID][chromosomeOffset].kMer[chromosomeItem[cip]].next.push_back(next);
			}
//			else
//				cout << "keep [" << chromosomeID << ", " << chromosomeOffset << ", " << cip << ", " << np << "] to [" << nextID << ", " << nextOffset << ", " << nip << "]" << endl;
		}
}

void switchSeqs(vector<Position> & p1, vector<Position> & p2, vector<char> & n1, vector<char> & n2)
{
	vector<Position> tp;
	vector<char> tn;

	tp = p1; p1 = p2; p2 = tp;
	tn = n1; n1 = n2; n2 = tn;
}

void updateGenomeWithRead(vector<Seq> & seqs, vector<vector<Base> > & genome, int length, int insertVariation, int mrl)
{
        int sp, pp, ppp, npp, gp, sp0, ip, ipp, np, tag0, tag1, tag2;// sp0 is for the other read
        unsigned int chromosomeID, chromosomeOffset, nextID, nextOffset, nextItem, cpp, i, chromosomeID0, chromosomeOffset0, nextID0, nextOffset0, size;
        Base b;
        ContiMer p;
	KMer k;
	vector<char> s, nextS;

        for(sp = 0, sp0 = 1; sp < seqs.size(); sp = sp + 2, sp0 = sp0 + 2)
	{
		pp = 0;
cont:
                for(; pp < seqs[sp].positionSets.size(); pp ++)
                {
                        for(ppp = 0; ppp < pp; ppp ++)
                                if(abs((int)(seqs[sp].positionSets[pp][0].chromosomeOffset - seqs[sp].positionSets[ppp][0].chromosomeOffset)) < seqs[sp].nucleotides.size())
				{
					pp ++;
                                        goto cont;
				}
			tag0 = tag1 = tag2 = 0;
                        if(seqs[sp].frSets[pp] == 1 && seqs[sp0].frSets[pp] == 0)
                        {
                                reverseComplement(seqs[sp].nucleotides);
				tag0 = 1;
			}
			else if(seqs[sp0].frSets[pp] == 1 && seqs[sp].frSets[pp] == 0)
			{
				reverseComplement(seqs[sp0].nucleotides);
				tag1 = 1;
                        }
			else
			{
				cout << "BOWTIE ALIGNMENT ERROR" << endl;
				exit(-1);
			}
			for(ppp = 0; ppp < seqs[sp].positionSets[pp].size() - length; ppp ++)
				if(seqs[sp].positionSets[pp][ppp].chromosomeID != -1 && seqs[sp0].positionSets[pp][ppp].chromosomeID != -1 &&
				seqs[sp].positionSets[pp][ppp].chromosomeOffset > seqs[sp0].positionSets[pp][ppp].chromosomeOffset)
				{
					switchSeqs(seqs[sp].positionSets[pp], seqs[sp0].positionSets[pp], seqs[sp].nucleotides, seqs[sp0].nucleotides);
					tag2 = 1;
					break;
				}
//has to always adjust nucleotides according to fr to accomodate reverse complements
                        for(ppp = 0; ppp < seqs[sp].positionSets[pp].size() - length; ppp ++)
                        {
                                if(seqs[sp].positionSets[pp][ppp].chromosomeID != -1)
                                {
                                        chromosomeID = seqs[sp].positionSets[pp][ppp].chromosomeID;
                                        chromosomeOffset = seqs[sp].positionSets[pp][ppp].chromosomeOffset;
					chromosomeID0 = seqs[sp0].positionSets[pp][ppp].chromosomeID;
					chromosomeOffset0 = seqs[sp0].positionSets[pp][ppp].chromosomeOffset;

                                        nextID = seqs[sp].positionSets[pp][ppp + 1].chromosomeID;
                                        nextOffset = seqs[sp].positionSets[pp][ppp + 1].chromosomeOffset;
					nextID0 = seqs[sp0].positionSets[pp][ppp + 1].chromosomeID;
					nextOffset0 = seqs[sp0].positionSets[pp][ppp + 1].chromosomeOffset;
//                                      nextItem = genome[nextID][nextOffset].nextPosition.size();
                                        if(nextID == -1)// insertion to genome
                                        {
                                                for(npp = ppp + 2; npp < seqs[sp].positionSets[pp].size(); npp ++)
                                                {
                                                        if(seqs[sp].positionSets[pp][npp].chromosomeID != -1)
                                                        {
                                                                nextID = seqs[sp].positionSets[pp][npp].chromosomeID;
                                                                nextOffset = seqs[sp].positionSets[pp][npp].chromosomeOffset;
								nextID0 = seqs[sp0].positionSets[pp][npp].chromosomeID;
								nextOffset0 = seqs[sp0].positionSets[pp][npp].chromosomeOffset;
                                                                if(seqs[sp].positionSets[pp][npp].chromosomeID == seqs[sp].positionSets[pp][ppp].chromosomeID && npp - ppp < MAX)// small insertion
                                                                {
									if(nextOffset == chromosomeOffset + 1)
//                                                                      for(cpp = chromosomeOffset; cpp < nextOffset; cpp ++)
                                                                        {
//                                                                                p.nextID = nextID;
//                                                                                p.nextOffset = nextOffset;
//                                                                                p.nextItem = genome[p.nextID][p.nextOffset].nextPosition.size();
//                                                                                p.contigOffset = ppp;
//                                                                                genome[chromosomeID][chromosomeOffset].nextPosition.push_back(p);
//										cout << "SI" << endl;
										if(s.size() > 0) s.clear();
										for(np = ppp; np < ppp + length; np ++)
											s.push_back(seqs[sp].nucleotides[np]);
										if(nextS.size() > 0) nextS.clear();
										//for(np = ppp + 1; np < ppp + length + 1; np ++)
										size = npp + length < seqs[sp].nucleotides.size() ? npp + length : seqs[sp].nucleotides.size();
										for(np = npp; np < size; np ++)
											nextS.push_back(seqs[sp].nucleotides[np]);
										updateKMer(genome, chromosomeID, chromosomeOffset, nextID, nextOffset, chromosomeID0, chromosomeOffset0, nextID0, nextOffset0, s, nextS, insertVariation, mrl);
                                                                                //genome[chromosomeID][cpp].nextID = chromosomeID;
                                                                                //genome[chromosomeID][cpp].nextOffset = cpp + 1;
                                                                        }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
									else
									{
			                                                        if(s.size() > 0) s.clear();
			                                                        for(np = ppp; np < ppp + length; np ++)
			                                                                s.push_back(seqs[sp].nucleotides[np]);
			                                                        if(nextS.size() > 0) nextS.clear();
			                                                        updateKMer(genome, chromosomeID, chromosomeOffset, chromosomeID, chromosomeOffset + 1, chromosomeID0, chromosomeOffset0, -1, -1, s, nextS, insertVariation, mrl);
			                                                        for(cpp = chromosomeOffset + 1; cpp < nextOffset - 1; cpp ++)
			                                                        {
			                                                                if(s.size() > 0) s.clear();
			                                                                if(nextS.size() > 0) nextS.clear();
			                                                                updateKMer(genome, chromosomeID, cpp, chromosomeID, cpp + 1, -1, -1, -1, -1, s, nextS, insertVariation, mrl);
			                                                        }
			                                                        if(s.size() > 0) s.clear();
			                                                        if(nextS.size() > 0) nextS.clear();
			                                                        //for(np = ppp + 1; np < ppp + length + 1; np ++)
										size = npp + length < seqs[sp].nucleotides.size() ? npp + length : seqs[sp].nucleotides.size();
										for(np = npp; np < size; np ++)
			                                                                nextS.push_back(seqs[sp].nucleotides[np]);
			                                                        updateKMer(genome, chromosomeID, cpp, chromosomeID, cpp + 1, -1, -1, nextID0, nextOffset0, s, nextS, insertVariation, mrl);

									}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                                                }
                                                                else// large insertion
//large insertion is ignored for simplicity
                                                                {
//                                                                      cout << "LI" << endl;
                                                                        p.nextID = chromosomeID;
                                                                        p.nextOffset = genome[chromosomeID].size();
                                                                        p.nextItem = 0;
                                                                        genome[chromosomeID][chromosomeOffset].contiMer.push_back(p);
                                                                        //genome[chromosomeID][chromosomeOffset].nextID = chromosomeID;
                                                                        //genome[chromosomeID][chromosomeOffset].nextOffset = genome[chromosomeID].size();
                                                                        for(i = 0; i < npp - ppp - 2; i ++)
                                                                        {
                                                                                b.contiMer.clear();
                                                                                b.nucleotide = seqs[sp].nucleotides[ppp + 1 + i];
                                                                                p.nextID = chromosomeID;
                                                                                p.nextOffset = genome[chromosomeID].size() + 1;
                                                                                p.nextItem = 0;
                                                                                b.contiMer.push_back(p);
                                                                                //b.nextID = chromosomeID;
                                                                                //b.nextOffset = genome[chromosomeID].size();
                                                                                genome[chromosomeID].push_back(b);
                                                                        }
                                                                        b.contiMer.clear();
                                                                        b.nucleotide = seqs[sp].nucleotides[npp];
                                                                        p.nextID = nextID; //contigs[cp].positionSets[pp][npp].chromosomeID;
                                                                        p.nextOffset = nextOffset; //contigs[cp].positionSets[pp][npp].chromosomeOffset;
                                                                        p.nextItem = genome[p.nextID][p.nextOffset].contiMer.size();
                                                                        b.contiMer.push_back(p);
                                                                        genome[chromosomeID].push_back(b);
                                                                        //genome[chromosomeID][genome[chromosomeID].size() - 1].nextID = nextID;
                                                                        //genome[chromosomeID][genome[chromosomeID].size() - 1].nextOffset = nextOffset;
                                                                }
                                                                ppp = npp - 1;
                                                                break;
                                                        }
                                                }
                                        }
                                        else if(nextID == chromosomeID && nextOffset != chromosomeOffset + 1)// deletion from genome
                                        {
                                                if(nextOffset - chromosomeOffset < SD)// small deletion
                                                {
//                                                      cout << "SD" << endl;
							if(s.size() > 0) s.clear();
							for(np = ppp; np < ppp + length; np ++)
								s.push_back(seqs[sp].nucleotides[np]);
							if(nextS.size() > 0) nextS.clear();
//							for(np = ppp + 1; np < ppp + length + 1; np ++)
//								nextS.push_back(seqs[sp].nucleotides[np]);
							updateKMer(genome, chromosomeID, chromosomeOffset, chromosomeID, chromosomeOffset + 1, chromosomeID0, chromosomeOffset0, -1, -1, s, nextS, insertVariation, mrl);
                                                        for(cpp = chromosomeOffset + 1; cpp < nextOffset - 1; cpp ++)
                                                        {
//                                                                p.nextID = chromosomeID;
//                                                                p.nextOffset = cpp + 1;
//                                                                p.nextItem = genome[p.nextID][p.nextOffset].nextPosition.size();
//                                                                genome[chromosomeID][cpp].nextPosition.push_back(p);
								if(s.size() > 0) s.clear();
								if(nextS.size() > 0) nextS.clear();
								updateKMer(genome, chromosomeID, cpp, chromosomeID, cpp + 1, -1, -1, -1, -1, s, nextS, insertVariation, mrl);
                                                                //genome[chromosomeID][cpp].nextID = chromosomeID;
                                                                //genome[chromosomeID][cpp].nextOffset = cpp + 1;
                                                        }
							if(s.size() > 0) s.clear();
							if(nextS.size() > 0) nextS.clear();
							for(np = ppp + 1; np < ppp + length + 1; np ++)
								nextS.push_back(seqs[sp].nucleotides[np]);
							updateKMer(genome, chromosomeID, cpp, chromosomeID, cpp + 1, -1, -1, nextID0, nextOffset0, s, nextS, insertVariation, mrl);
//The empty s and nextS make it possible to generate single nucleotide contig (the min contig size is not k-mer size)
                                                }
                                                else// large deletion
                                                {
//                                                        cout << "LD" << endl;
//                                                        p.nextID = nextID;
//                                                        p.nextOffset = nextOffset;
//                                                        p.nextItem = genome[p.nextID][p.nextOffset].nextPosition.size();
//                                                        genome[chromosomeID][chromosomeOffset].nextPosition.push_back(p);
							if(s.size() > 0) s.clear();
							for(np = ppp; np < ppp + length; np ++)
								s.push_back(seqs[sp].nucleotides[np]);
							if(nextS.size() > 0) nextS.clear();
							for(np = ppp + 1; np < ppp + length + 1; np ++)
								nextS.push_back(seqs[sp].nucleotides[np]);
							updateKMer(genome, chromosomeID, chromosomeOffset, nextID, nextOffset, chromosomeID0, chromosomeOffset0, nextID0, nextOffset0, s, nextS, insertVariation, mrl);
                                                        //genome[chromosomeID][chromosomeOffset].nextID = nextID;
                                                        //genome[chromosomeID][chromosomeOffset].nextOffset = nextOffset;

                                                }
                                        }
                                        else// ordinary case
                                        {
//						  cout << "OD" << endl;
//                                                p.nextID = nextID;
//                                                p.nextOffset = nextOffset;
//                                                p.nextItem = genome[p.nextID][p.nextOffset].nextPosition.size();
//                                                genome[chromosomeID][chromosomeOffset].nextPosition.push_back(p);
						if(s.size() > 0) s.clear();
						for(np = ppp; np < ppp + length; np ++)
							s.push_back(seqs[sp].nucleotides[np]);
						if(nextS.size() > 0) nextS.clear();
						for(np = ppp + 1; np < ppp + length + 1; np ++)
							nextS.push_back(seqs[sp].nucleotides[np]);
						updateKMer(genome, chromosomeID, chromosomeOffset, nextID, nextOffset, chromosomeID0, chromosomeOffset0, nextID0, nextOffset0, s, nextS, insertVariation, mrl);
                                                //genome[chromosomeID][chromosomeOffset].nextID = nextID;
                                                //genome[chromosomeID][chromosomeOffset].nextOffset = nextOffset;
                                        }
                                }
                        }
			if(tag2 == 1) switchSeqs(seqs[sp].positionSets[pp], seqs[sp0].positionSets[pp], seqs[sp].nucleotides, seqs[sp0].nucleotides);
			if(tag0 == 1 && tag1 == 0) reverseComplement(seqs[sp].nucleotides);
			else if(tag1 == 1 && tag0 == 0) reverseComplement(seqs[sp0].nucleotides);
			else
			{
				cout << "UNKNOWN ERROR" << endl;
				exit(-1);
			}
///                        if(nextID != -1)
///                        {
//                                p.nextID = p.nextOffset = p.nextItem = -1;
//                                genome[nextID][nextOffset].nextPosition.push_back(p);
///				updateKMer(genome, nextID, nextOffset, -1, -1, nextID0, nextOffset0, -1, -1);
///                        }
///                        else// chromosomeID != -1 and nextID == -1
///                        {
//                                p.nextID = p.nextOffset = p.nextItem = -1;
//                                genome[chromosomeID][chromosomeOffset].nextPosition.push_back(p);
///				updateKMer(genome, chromosomeID, chromosomeOffset, -1, -1, chromosomeID0, chromosomeOffset0, -1, -1);
///                        }
                }
	}
}

void loadReadAlignment(vector<vector<Base> > & genome, int length, int insertVariation, int chromosomeID, int mrl)
{ 
	vector<Seq> reads;
	ifstream r, ra;
	int seqStartID, aliStartID, finish, rp;
	string s;
//	int counter = 0;

	r.open("tmp/_reads.fa");
	s = "tmp/_reads_genome." + itoa(chromosomeID) + ".bowtie";
	ra.open(s.c_str());
	
	seqStartID = aliStartID = -1;
cont:
        finish = loadSeq(r, reads, aliStartID, seqStartID);
        loadReadAli(ra, reads, aliStartID, seqStartID);
        updateGenomeWithRead(reads, genome, length, insertVariation, mrl);
	if(finish == 0)// && counter < 2)
	{
//		cout << counter ++ << endl;
		reads.clear();
		goto cont;
	}
}

typedef struct contigStruct
{
	int extended;
///	unsigned int contigID;
///	unsigned int contigID0;
///	unsigned int contigOffset;
///	unsigned int contigOffset0;
	unsigned int startID;
	unsigned int startOffset;
	unsigned int endID;
	unsigned int endOffset;
	unsigned int startID0;
	unsigned int startOffset0;
	unsigned int endID0;
	unsigned int endOffset0;
	vector<char> nucleotides;
	/*vector<char> pairedNucleotides;*/
} Contig;

vector<Contig> contigs;

int contain(unsigned int startID1, unsigned int startOffset1, unsigned int endID1, unsigned int endOffset1, unsigned int startID2, unsigned int startOffset2, unsigned int endID2, unsigned int endOffset2)
{
	if(startID1 == startID2 && endID1 == endID2 && startOffset1 <= startOffset2 && endOffset1 >= endOffset2)
		return 1;
	return 0;
}

void filterLowCoverage(vector<vector<Base> > & genome, int coverage)
{
	unsigned int gp, cp, ip, cpp, pp;

//check point
/*
          cout << "*************************" << endl;
          for(gp = 0; gp < genome.size(); gp ++)
          {
                  for(cpp = 0; cpp < genome[gp].size(); cpp ++)
                  {
                          cout << "[";
                          for(pp = 0; pp < genome[gp][cpp].contiMer.size(); pp ++)
                                  cout << "<" << genome[gp][cpp].contiMer[pp].contigID << ", " << genome[gp][cpp].contiMer[pp].contigOffset << "| " << genome[gp][cpp].contiMer[pp].previousID << ", " << genome[gp][cpp].contiMer[pp].previousOffset << ", " << genome[gp][cpp].contiMer[pp].previousItem << "| " << genome[gp][cpp].contiMer[pp].nextID << ", " << genome[gp][cpp].contiMer[pp].nextOffset << ", " << genome[gp][cpp].contiMer[pp].nextItem << ">";
                          cout << "]";
                        if((cpp + 1) % 60 == 0) cout << endl;
                  }
                  cout << endl;
          }
        cout << "-------------------------" << endl;

        Next next;
        Previous previous;
        next.nextID = next.nextOffset = next.nextItem = -1;
        previous.previousID = previous.previousOffset = previous.previousItem = -1;
        for(gp = 0; gp < genome.size(); gp ++)
                for(cpp = 0; cpp < genome[gp].size(); cpp ++)
                        for(ip = 0; ip < genome[gp][cpp].kMer.size(); ip ++)
                        {
                                if(genome[gp][cpp].kMer[ip].next.size() == 0)
                                        genome[gp][cpp].kMer[ip].next.push_back(next);
                                //if(genome[gp][cpp].kMer[ip].previous.size() == 0)
                                        //genome[gp][cpp].kMer[ip].previous.push_back(previous);
                        }

          for(gp = 0; gp < genome.size(); gp ++)
          {
                  for(cpp = 0; cpp < genome[gp].size(); cpp ++)
                  {
                          cout << "[";
                          for(pp = 0; pp < genome[gp][cpp].kMer.size(); pp ++)
				if(genome[gp][cpp].kMer[pp].traversed == 0)
                                	cout << "<" << genome[gp][cpp].kMer[pp].contigID << ", " << genome[gp][cpp].kMer[pp].contigOffset << ", " << genome[gp][cpp].kMer[pp].contigID0 << ", " << genome[gp][cpp].kMer[pp].contigOffset0 << ", " << genome[gp][cpp].kMer[pp].chromosomeID0 << ", " << genome[gp][cpp].kMer[pp].chromosomeOffset0 <<  "| " << genome[gp][cpp].kMer[pp].previous[0].previousID << ", " << genome[gp][cpp].kMer[pp].previous[0].previousOffset << ", " << genome[gp][cpp].kMer[pp].previous[0].previousItem << "| " << genome[gp][cpp].kMer[pp].next[0].nextID << ", " << genome[gp][cpp].kMer[pp].next[0].nextOffset << ", " << genome[gp][cpp].kMer[pp].next[0].nextItem << ">";
                          cout << "]";
                        if((cpp + 1) % 60 == 0) cout << endl;
                  }
          }
          cout << endl;
*/
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	for(gp = 0; gp < genome.size(); gp ++)
		for(cp = 0; cp < genome[gp].size(); cp ++)
			for(ip = 0; ip < genome[gp][cp].kMer.size(); ip ++)
			{
				if(genome[gp][cp].kMer[ip].contigID == -1)
				{
					if(genome[gp][cp].kMer[ip].coverage < coverage)
						genome[gp][cp].kMer[ip].traversed = 1;
				}
//				else
//				{
//					if(genome[gp][cp].kMer[ip].coverage < coverage / 5)
//						genome[gp][cp].kMer[ip].traversed = 1;
//				}
			}
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/*
        cout << "-------------------------" << endl;

          for(gp = 0; gp < genome.size(); gp ++)
          {
                  for(cpp = 0; cpp < genome[gp].size(); cpp ++)
                  {
                          cout << "[";
                          for(pp = 0; pp < genome[gp][cpp].kMer.size(); pp ++)
				if(genome[gp][cpp].kMer[pp].traversed == 0)
                                	cout << "<" << genome[gp][cpp].kMer[pp].contigID << ", " << genome[gp][cpp].kMer[pp].contigOffset << ", " << genome[gp][cpp].kMer[pp].contigID0 << ", " << genome[gp][cpp].kMer[pp].contigOffset0 << ", " << genome[gp][cpp].kMer[pp].chromosomeID0 << ", " << genome[gp][cpp].kMer[pp].chromosomeOffset0 <<  "| " << genome[gp][cpp].kMer[pp].previous[0].previousID << ", " << genome[gp][cpp].kMer[pp].previous[0].previousOffset << ", " << genome[gp][cpp].kMer[pp].previous[0].previousItem << "| " << genome[gp][cpp].kMer[pp].next[0].nextID << ", " << genome[gp][cpp].kMer[pp].next[0].nextOffset << ", " << genome[gp][cpp].kMer[pp].next[0].nextItem << ">";
                          cout << "]";
                        if((cpp + 1) % 60 == 0) cout << endl;
                  }
          }
          cout << endl;
          cout << "*************************" << endl;
*/

}

void countBranches(vector<vector<Base> > & genome)
{
	unsigned int gp, cp, ip, nBranches = 0;

	for(gp = 0; gp < genome.size(); gp ++)
		for(cp = 0; cp < genome[gp].size(); cp ++)
			for(ip = 0; ip < genome[gp][cp].kMer.size(); ip ++)
				if(genome[gp][cp].kMer[ip].traversed == 0 && genome[gp][cp].kMer[ip].next.size() > 1)
					nBranches ++;

	cout << nBranches << " branches in de Bruijn graph" << endl;
}

void checkContiMers(vector<vector<Base> > & genome)
{
        unsigned int gp, cp, ip, nBranches = 0;

        for(gp = 0; gp < genome.size(); gp ++)
                for(cp = 0; cp < genome[gp].size(); cp ++)
                        for(ip = 0; ip < genome[gp][cp].contiMer.size(); ip ++)
                                if(genome[gp][cp].contiMer[ip].contigID != -1 && genome[gp][cp].contiMer[ip].contigOffset == -1)
                                        cout << "error @" << gp << ", " << cp << ", " << ip << endl;

}

char max(int a, int c, int g, int t, int n)
{
	if(a == 0 && c == 0 && g == 0 && t == 0 && n == 0) return 'X';
	if(a >= c && a >= g && a >= t && a >= n) return 'A';
	if(c >= a && c >= g && c >= t && c >= n) return 'C';
	if(g >= a && g >= c && g >= t && g >= n) return 'G';
	if(t >= a && t >= c && t >= g && t >= n) return 'T';
	if(n >= a && n >= c && n >= g && n >= t) return 'N';
}

//string getKMer(vector<char> A, vector<char> C, vector<char> G, vector<char> T, vector<char> N)
//{
//	int np;
//	string s;
//
//	for(np = 0; np < A.size(); np ++)
//		s.push_back(max(A[np], C[np], G[np], T[np], N[np]));
//	return s;
//}

void extdContigs1(vector<vector<Base> > & genome, int coverage, int k, int chromosomeID)
{
	unsigned int gp, cp, ip, gpp, cpp, ipp, gppBak, cppBak, ippBak, seqID, i, kMerTag, numItem, nextItem, sp, gpBak, cpBak, ipBak, np, gpp0, cpp0, ngp, ncp, nip, ippp, item, next, count, nCount, startIDBak = -1, startOffsetBak = -1, endIDBak = -1, endOffsetBak = -1;
	Contig contig;
//	vector<Contig> contigs;
	vector<char> sBak;
	vector<Contig> contigsBuf;
	char nucleotide;

//******************************************************************************************************************************************************************************
	ofstream out;
        string s = "tmp/_pre_extended_contigs." + itoa(chromosomeID) + ".fa";
        out.open(s.c_str());
//******************************************************************************************************************************************************************************

	filterLowCoverage(genome, coverage);
//	countBranches(genome);

	seqID = 0;
	for(gp = 0; gp < genome.size(); gp ++)
		for(cp = 0; cp < genome[gp].size();)
		{
//			if(genome[gp][cp].kMer.size() != 0)
				for(ip = 0; ip < genome[gp][cp].kMer.size(); ip ++)
				{
					if(genome[gp][cp].kMer[ip].traversed == 0)
					{
//						contig.clear();
						gpp = gp;
						cpp = cp;
						ipp = ip;
						kMerTag = 1;
						contigs.push_back(contig);
						contigs[contigs.size() - 1].startID = gp;
						contigs[contigs.size() - 1].startOffset = cp;
						contigs[contigs.size() - 1].startID0 = genome[gp][cp].kMer[ip].chromosomeID0;
						contigs[contigs.size() - 1].startOffset0 = genome[gp][cp].kMer[ip].chromosomeOffset0;
						contigs[contigs.size() - 1].extended = 0;//added
///						contigs[contigs.size() - 1].contigOffset = -1;

						while(kMerTag == 1 && genome[gpp][cpp].kMer[ipp].traversed == 0 || kMerTag == 0)
						{
							if(kMerTag == 0)
								contigs[contigs.size() - 1].nucleotides.push_back(genome[gpp][cpp].contiMer[ipp].nucleotide);
							else
							{
//								nucleotide = max(genome[gpp][cpp].kMer[ipp].A[0], genome[gpp][cpp].kMer[ipp].C[0], genome[gpp][cpp].kMer[ipp].G[0], genome[gpp][cpp].kMer[ipp].T[0], genome[gpp][cpp].kMer[ipp].N[0]);
								nucleotide = max(genome[gpp][cpp].kMer[ipp].A, genome[gpp][cpp].kMer[ipp].C, genome[gpp][cpp].kMer[ipp].G, genome[gpp][cpp].kMer[ipp].T, genome[gpp][cpp].kMer[ipp].N);
								if(nucleotide != 'X')
									contigs[contigs.size() - 1].nucleotides.push_back(nucleotide);
								else
									contigs[contigs.size() - 1].nucleotides.push_back(genome[gpp][cpp].nucleotide);
							}

							if(kMerTag == 1 && genome[gpp][cpp].kMer[ipp].contigOffset != -1 || kMerTag == 0) 
							//make the change: incorporate novel ones
								contigs[contigs.size() - 1].extended = 1;//added

							if(kMerTag == 1)
							{
								gpp0 = genome[gpp][cpp].kMer[ipp].chromosomeID0;
								cpp0 = genome[gpp][cpp].kMer[ipp].chromosomeOffset0;
								/*if(gpp0 != -1)
									contigs[contigs.size() - 1].pairedNucleotides.push_back(genome[gpp0][cpp0].nucleotide);*/

///								if(contigs[contigs.size() - 1].contigOffset == -1 && genome[gpp][cpp].kMer[ipp].contigOffset != -1)
///								{
///									contigs[contigs.size() - 1].contigID = genome[gpp][cpp].kMer[ipp].contigID;
///									contigs[contigs.size() - 1].contigID0 = genome[gpp][cpp].kMer[ipp].contigID0;
///									contigs[contigs.size() - 1].contigOffset = genome[gpp][cpp].kMer[ipp].contigOffset;
///									contigs[contigs.size() - 1].contigOffset0 = genome[gpp][cpp].kMer[ipp].contigOffset0;
///								}

//								cout << "k-mer " << cpp << ": "  << genome[gpp][cpp].nucleotide << " | ";
								genome[gpp][cpp].kMer[ipp].traversed = 1;
								gpBak = gpp;
								cpBak = cpp;
								ipBak = ipp;
								sBak = genome[gpp][cpp].kMer[ipp].s;//getKMer(genome[gpp][cpp].kMer[ipp].A, genome[gpp][cpp].kMer[ipp].C, genome[gpp][cpp].kMer[ipp].G, genome[gpp][cpp].kMer[ipp].T, genome[gpp][cpp].kMer[ipp].N);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

								nCount = 0;
								next = -1;
								for(np = 0; np < genome[gpp][cpp].kMer[ipp].next.size(); np ++)
								{
									ngp = genome[gpp][cpp].kMer[ipp].next[np].nextID;
									ncp = genome[gpp][cpp].kMer[ipp].next[np].nextOffset;
									nip = genome[gpp][cpp].kMer[ipp].next[np].nextItem;
									if(ngp != -1 && genome[ngp][ncp].kMer[nip].traversed == 0)
									{
										next = np;
										nCount ++;
									}
								}
								if(nCount == 1)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//								if(genome[gpp][cpp].kMer[ipp].next.size() == 1 && genome[gpp][cpp].kMer[ipp].next[0].nextID != -1)// need improvement considering filteration
								{
									gppBak = genome[gpp][cpp].kMer[ipp].next[next].nextID;
									cppBak = genome[gpp][cpp].kMer[ipp].next[next].nextOffset;
									ippBak = genome[gpp][cpp].kMer[ipp].next[next].nextItem;
									gpp = gppBak;
									cpp = cppBak;
									ipp = ippBak;
									kMerTag = 1;
//									cout << "keep k-mer" << " | ";
								}
								else if(genome[gpp][cpp].contiMer.size() == 1 && genome[gpp][cpp].contiMer[0].nextID != -1)
								{
									gppBak = genome[gpp][cpp].contiMer[0].nextID;
									cppBak = genome[gpp][cpp].contiMer[0].nextOffset;
									ippBak = genome[gpp][cpp].contiMer[0].nextItem;
									gpp = gppBak;
									cpp = cppBak;
									ipp = ippBak;
									kMerTag = 0;
//									cout << "switch from k-mer to conti-mer" << " | ";
								}
								else
									kMerTag = -1;
							}
							else// if(kMerTag == 0)
							{
///								if(contigs[contigs.size() - 1].contigOffset == -1)
///								{
///									contigs[contigs.size() - 1].contigOffset = genome[gpp][cpp].contiMer[ipp].contigOffset;
///									contigs[contigs.size() - 1].contigOffset0 = -1;
///								}
//								cout << "conti-mer " << cpp << ": " << genome[gpp][cpp].nucleotide << " | ";
								if(genome[gpp][cpp].contiMer[ipp].nextID != -1)
								{
									gppBak = genome[gpp][cpp].contiMer[ipp].nextID;
									cppBak = genome[gpp][cpp].contiMer[ipp].nextOffset;
									ippBak = genome[gpp][cpp].contiMer[ipp].nextItem;
									gpp = gppBak;
									cpp = cppBak;
									ipp = ippBak;
									kMerTag = 0;
//maybe not necessary to get back to k-mer?
//									for(ippBak = 0, numItem = 0; ippBak < genome[gpp][cpp].kMer.size(); ippBak ++)
//										if(genome[gpp][cpp].kMer[ippBak].traversed == 0)
//										{
//											numItem ++;
//											nextItem = ippBak;
//										}
//									if(numItem == 1)
//									{
//										cout << "switch from conti-mer to k-mer (1)" << " | ";
//										ipp = nextItem;
//										kMerTag = 1;
//									}
//									else
//										cout << "keep conti-mer" << " | ";
								}
								else
								{
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

									count = nCount = 0;
									item = next = -1;
									for(ippp = 0; ippp < genome[gpp][cpp].kMer.size(); ippp ++)
										if(genome[gpp][cpp].kMer[ippp].traversed == 0)
										{
											count ++;
											item = ippp;
										}
									if(count == 1)
									{
										for(np = 0; np < genome[gpp][cpp].kMer[item].next.size(); np ++)
										{
											ngp = genome[gpp][cpp].kMer[item].next[np].nextID;
											ncp = genome[gpp][cpp].kMer[item].next[np].nextOffset;
											nip = genome[gpp][cpp].kMer[item].next[np].nextItem;
											if(ngp != -1 && genome[ngp][ncp].kMer[nip].traversed == 0)
											{
												nCount ++;
												next = np;
											}
										}
									}
                                                                	if(nCount == 1)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//									if(genome[gpp][cpp].kMer.size() == 1 && genome[gpp][cpp].kMer[0].next.size() == 1 && genome[gpp][cpp].kMer[0].next[0].nextID != -1)// need improvement considering filteration
									{
										gppBak = genome[gpp][cpp].kMer[item].next[next].nextID;
										cppBak = genome[gpp][cpp].kMer[item].next[next].nextOffset;
										ippBak = genome[gpp][cpp].kMer[item].next[next].nextItem;
										gpp = gppBak;
										cpp = cppBak;
										ipp = ippBak;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										if(genome[gpp][cpp].kMer[ipp].traversed == 0)
											kMerTag = 1;
										else
											kMerTag = -2;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										//kMerTag = 1;
//										cout << "switch from conti-mer to k-mer (2)" << " | ";
									}
									else
										kMerTag = -2;
								}
							}
						}

//xpp and xppBak are always the same
						if(kMerTag == 1)// && genome[gpp][cpp].kMer[ipp].traversed == 1
						{
							contigs[contigs.size() - 1].endID = gppBak;
							contigs[contigs.size() - 1].endOffset = cppBak;
						}
						else//kMerTag == -1 || kMerTag == -2
						{
                                                	contigs[contigs.size() - 1].endID = gpp;
                                                	contigs[contigs.size() - 1].endOffset = cpp;
						}

//						if(kMerTag == 1)
//						{
//							contigs[contigs.size() - 1].endID0 = genome[gppBak][cppBak].kMer[ippBak].chromosomeID0;
//							contigs[contigs.size() - 1].endOffset0 = genome[gppBak][cppBak].kMer[ippBak].chromosomeOffset0;
//						}
//						else if(kMerTag == -1)
						if(kMerTag == 1 || kMerTag == -1)
						{
                                                	contigs[contigs.size() - 1].endID0 = genome[gpp][cpp].kMer[ipp].chromosomeID0;
                                                	contigs[contigs.size() - 1].endOffset0 = genome[gpp][cpp].kMer[ipp].chromosomeOffset0;
						}
						else// kMerTag == -2
						{
							contigs[contigs.size() - 1].endID0 = -1;
							contigs[contigs.size() - 1].endOffset0 = -1;
						}

//						cout << endl;
//						if(contigs[contigs.size() - 1].nucleotides.size() == 1)
//							cout << "[" << gpp << ", " << cpp << ", " << ipp << "]: " << genome[gpp][cpp].kMer[ipp].next.size() << endl;
						if(kMerTag == -1 || kMerTag == 1)
						{
//append the last |s| - 1 bases only if getting out of while from k-mer rather than conti-mer
							for(sp = 1; sp < sBak.size(); sp ++)
								contigs[contigs.size() - 1].nucleotides.push_back(sBak[sp]);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
							contigs[contigs.size() - 1].endOffset = contigs[contigs.size() - 1].endOffset + sBak.size() - 1;
							contigs[contigs.size() - 1].endOffset0 = contigs[contigs.size() - 1].endOffset0 + sBak.size() - 1;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						}
//						cout << endl;

//						cout << ">" << seqID ++ << endl;
//						for(i = 0; i < contig.size(); i ++)
//						{
//							cout << contig[i];
//							if((i + 1) % 60 == 0 || i == contig.size() - 1)
//								cout << endl;
//						}

//						if(contigs[contigs.size() - 1].extended == 1)
//						{
//							cout << ">" << seqID ++ << endl;
//							for(i = 0; i < contigs[contigs.size() - 1].nucleotides.size(); i ++)
//							{
//								cout << contigs[contigs.size() - 1].nucleotides[i];
//								if((i + 1) % 60 == 0 || i == contigs[contigs.size() - 1].nucleotides.size() - 1)
//									cout << endl;
//							}
//						}
//******************************************************************************************************************************************************************************
						if(/*contigs[contigs.size() - 1].extended == 1 && */ contain(startIDBak, startOffsetBak, endIDBak, endOffsetBak, contigs[contigs.size() - 1].startID, contigs[contigs.size() - 1].startOffset, contigs[contigs.size() - 1].endID, contigs[contigs.size() - 1].endOffset) == 0)
						{
							out << ">" << seqID ++ << ", " /*<< contigs[contigs.size() - 1].contigID << ", " << contigs[contigs.size() - 1].contigOffset << ", " << contigs[contigs.size() - 1].contigID0 << ", " << contigs[contigs.size() - 1].contigOffset0 << ", "*/ << contigs[contigs.size() - 1].extended << ", " << contigs[contigs.size() - 1].startID << ", " << contigs[contigs.size() - 1].startOffset << ", " << contigs[contigs.size() - 1].endID << ", " << contigs[contigs.size() - 1].endOffset << ", " << contigs[contigs.size() - 1].startID0 << ", " << contigs[contigs.size() - 1].startOffset0 << ", " << contigs[contigs.size() - 1].endID0 << ", " << contigs[contigs.size() - 1].endOffset0 << " " << endl;
							for(i = 0; i < contigs[contigs.size() - 1].nucleotides.size(); i ++)
							{
								out << contigs[contigs.size() - 1].nucleotides[i];
								if((i + 1) % 60 == 0 || i == contigs[contigs.size() - 1].nucleotides.size() - 1)
									out << endl;
							}
							startIDBak = contigs[contigs.size() - 1].startID;
							startOffsetBak = contigs[contigs.size() - 1].startOffset;
							endIDBak = contigs[contigs.size() - 1].endID;
							endOffsetBak = contigs[contigs.size() - 1].endOffset;
						}
						contigs.clear();
//******************************************************************************************************************************************************************************
					}
				}
			if(endOffsetBak - startOffsetBak > 100000)
			{
				if(gp == endIDBak && cp + 1000 < endOffsetBak)
					cp = cp + 1000;
				else
					cp ++;
			}
			else
				cp ++;
		}
}

void parse(string buf, int & extended, /*unsigned int & contigID, unsigned int & contigOffset, unsigned int & contigID0, unsigned int & contigOffset0,*/  unsigned int & startID, unsigned int & startOffset, unsigned int & endID, unsigned int & endOffset, unsigned int & startID0, unsigned int & startOffset0, unsigned int & endID0, unsigned int & endOffset0)
{
	int i, item = 0, j1 = 0, j2 = 0, j3 = 0, j4 = 0, j5 = 0, j6 = 0, j7 = 0, j8 = 0, j9 = 0, j10 = 0, j11 = 0, j12 = 0, j13 = 0;
	char contigIDBuf[20] = {'\n'}, contigOffsetBuf[20] = {'\n'}, contigID0Buf[20] = {'\n'}, contigOffset0Buf[20] = {'\n'}, extendedBuf[20] = {'\n'}, startIDBuf[20] = {'\n'}, startOffsetBuf[20] = {'\n'}, endIDBuf[20] = {'\n'}, endOffsetBuf[20] = {'\n'}, startID0Buf[20] = {'\n'}, startOffset0Buf[20] = {'\n'}, endID0Buf[20] = {'\n'}, endOffset0Buf[20] = {'\n'};

	for(i = 0; i < buf.size(); i ++)
	{
		if(buf[i] == ' ')
		{
			item ++;
			continue;
		}

		if(item == 0)
			continue;
///		else if(item == 1)
///			contigIDBuf[j1 ++] = buf[i];
///		else if(item == 2)
///			contigOffsetBuf[j2 ++] = buf[i];
///		else if(item == 3)
///			contigID0Buf[j3 ++] = buf[i];
///		else if(item == 4)
///			contigOffset0Buf[j4 ++] = buf[i];
		else if(item == 1)
			extendedBuf[j1 ++] = buf[i];
		else if(item == 2)
			startIDBuf[j2 ++] = buf[i];
		else if(item == 3)
			startOffsetBuf[j3 ++] = buf[i];
		else if(item == 4)
			endIDBuf[j4 ++] = buf[i];
		else if(item == 5)
			endOffsetBuf[j5 ++] = buf[i];
                else if(item == 6)
                        startID0Buf[j6 ++] = buf[i];
                else if(item == 7)
                        startOffset0Buf[j7 ++] = buf[i];
                else if(item == 8)
                        endID0Buf[j8 ++] = buf[i];
                else if(item == 9)
                        endOffset0Buf[j9 ++] = buf[i];
		else if(item == 10)
			break;
		else
		{
			cout << "UNKNOWN ERROR" << endl;
			exit(-1);
		}
	}
	extended = atoi(extendedBuf);
///	contigID = atoi(contigIDBuf);
///	contigOffset = atoi(contigOffsetBuf);
///	contigID0 = atoi(contigID0Buf);
///	contigOffset0 = atoi(contigOffset0Buf);
	startID = atoi(startIDBuf);
	startOffset = atoi(startOffsetBuf);
	endID = atoi(endIDBuf);
	endOffset = atoi(endOffsetBuf);
        startID0 = atoi(startID0Buf);
        startOffset0 = atoi(startOffset0Buf);
        endID0 = atoi(endID0Buf);
        endOffset0 = atoi(endOffset0Buf);
}

void loadContigs(int chromosomeID)
{
	int cp, cpp, i;
	ifstream in;
	string buf;
	Contig contig;

	string s = "tmp/_pre_extended_contigs." + itoa(chromosomeID) + ".fa";
        in.open(s.c_str());
	contig.extended = /*contig.contigID = contig.contigOffset = contig.contigID0 = contig.contigOffset0 =*/ contig.startID = contig.startOffset = contig.endID = contig.endOffset = contig.startID0 = contig.startOffset0 = contig.endID0 = contig.endOffset0 = -1;

	if(in.is_open())
	{
		while(in.good())
		{
			getline(in, buf);
			if(buf[0] == 0)
				break;

			if(buf[0] == '>')
			{
				contigs.push_back(contig);
				parse(buf, contigs[contigs.size() - 1].extended, /*contigs[contigs.size() - 1].contigID, contigs[contigs.size() - 1].contigOffset, contigs[contigs.size() - 1].contigID0, contigs[contigs.size() - 1].contigOffset0,*/ contigs[contigs.size() - 1].startID, contigs[contigs.size() - 1].startOffset, contigs[contigs.size() - 1].endID, contigs[contigs.size() - 1].endOffset, contigs[contigs.size() - 1].startID0, contigs[contigs.size() - 1].startOffset0, contigs[contigs.size() - 1].endID0, contigs[contigs.size() - 1].endOffset0);
			}
			else
			{
				for(i = 0; i < buf.size(); i ++)
					contigs[contigs.size() - 1].nucleotides.push_back(buf[i]);
			}
		}
	}
	else
	{
		cout << "CANNOT OPEN FILE!" << endl;
		exit(-1);
	}
}

void extdContigs2(int chromosomeID)
{
	int cp, cpp, cppBak, np;
	vector<Contig> contigsBuf;

	loadContigs(chromosomeID);

	for(cp = 0; cp < contigs.size(); cp ++)
	{
		if(contigs[cp].extended == 1)
		{
			for(cpp = cp + 1; cpp < contigs.size(); cpp ++)
			//if(contigs[cpp].extended == 1)
			{
				if(contain(contigs[cp].startID, contigs[cp].startOffset, contigs[cp].endID, contigs[cp].endOffset,
				contigs[cpp].startID, contigs[cpp].startOffset, contigs[cpp].endID, contigs[cpp].endOffset))
				{
//					cout << "duplication removed" << endl;
					contigs[cpp].extended = 2;
				}
				else if(contigs[cp].endID != contigs[cpp].startID || contigs[cp].endOffset < contigs[cpp].startOffset)
					break;
			}
		}
	}

	for(cp = contigs.size() - 1; cp != -1; cp --)
	{
		if(contigs[cp].extended == 1)
		{
			for(cpp = cp - 1; cpp != -1; cpp --)
			//if(contigs[cpp].extended == 1)
			{
				if(contain(contigs[cp].startID, contigs[cp].startOffset, contigs[cp].endID, contigs[cp].endOffset,
                                contigs[cpp].startID, contigs[cpp].startOffset, contigs[cpp].endID, contigs[cpp].endOffset))
				{
//					cout << "duplication removed" << endl;
					contigs[cpp].extended = 2;
				}
				else if(contigs[cpp].endID != contigs[cp].startID || contigs[cpp].endOffset < contigs[cp].startOffset)
					break;
			}
		}
	}
//remove duplication (0: no contig; 1: extended with contig; 2: filtered)

	for(cp = 0; cp < contigs.size(); cp ++)
	{
cont:
		if(contigs[cp].extended == 1)
		{
			cppBak = -1;
			contigsBuf.clear();
			for(cpp = cp + 1; cpp < contigs.size(); cpp ++)
				if(contigs[cpp].extended != 2)
				{
					if(//contigs[cp].endID == contigs[cpp].startID && 
					//contigs[cp].endOffset < contigs[cpp].startOffset && 
					contigs[cp].endOffset >= contigs[cpp].startOffset)//
					{
//						cout << "potential extension: " << contigs[cp].startOffset << " to " << contigs[cpp].startOffset << endl;
						contigsBuf.push_back(contigs[cpp]);
						cppBak = cpp;
					}
					else if(//contigs[cp].endID != contigs[cpp].startID || 
					contigs[cp].endOffset < contigs[cpp].startOffset)
						break;
				}
			if(contigsBuf.size() == 1)
			{
//				cout << "real extension: " << contigs[cp].startOffset << " to " <<  contigs[cppBak].startOffset << endl;
				contigs[cppBak].extended = 2; 
				for(np = //contigs[cp].endOffset < contigs[cpp].startOffset ? 0 : //
					contigs[cp].endOffset - contigsBuf[0].startOffset + 1; np < contigsBuf[0].nucleotides.size(); np ++)
					contigs[cp].nucleotides.push_back(contigsBuf[0].nucleotides[np]);
				/*for(np = //contigs[cp].endOffset < contigs[cpp].startOffset ? 0 : //
					contigs[cp].endOffset - contigs[cpp].startOffset + 1; np < contigsBuf[0].pairedNucleotides.size(); np ++)
					contigs[cp].pairedNucleotides.push_back(contigsBuf[0].pairedNucleotides[np]);*/
				contigs[cp].endID = contigsBuf[0].endID;
                                contigs[cp].endOffset = contigsBuf[0].endOffset;
                                contigs[cp].endID0 = contigsBuf[0].endID0;
                                contigs[cp].endOffset0 = contigsBuf[0].endOffset0;
				goto cont;
			}
		}
	}
//join contigs

/*
	for(cp = 0; cp < contigs.size(); cp ++)
	{
		if(contigs[cp].extended == 1 && contigs[cp].startOffset0 != -1)
		// && contigs[cp].endID == contigs[cp].startID0)// not necessary for single genome loading
		{
			if(contigs[cp].nucleotides.size() < contigs[cp].startOffset0 - contigs[cp].startOffset)
			{
				continue;
//				for(np = 0; np < contigs[cp].startOffset0 - contigs[cp].startOffset - contigs[cp].nucleotides.size(); np ++)
//					contigs[cp].nucleotides.push_back('N');// possible for single genome loading
				for(np = contigs[cp].startOffset + contigs[cp].nucleotides.size(); np < contigs[cp].startOffset0; np ++)
					contigs[cp].nucleotides.push_back(genome[0][np].nucleotide);
				for(np = 0; np < contigs[cp].pairedNucleotides.size(); np ++)
					contigs[cp].nucleotides.push_back(contigs[cp].pairedNucleotides[np]);
			}
			else
			{
				for(np = contigs[cp].nucleotides.size() - (contigs[cp].startOffset0 - contigs[cp].startOffset); np < contigs[cp].pairedNucleotides.size(); np ++)
					contigs[cp].nucleotides.push_back(contigs[cp].pairedNucleotides[np]);
			}
		}
	}
*/
//scaffolding

//	cout << endl;
//	int seqID = 0;
//	for(cp = 0; cp < contigs.size(); cp ++)
//	{
//		if(contigs[cp].extended == 1)
//		{
//			cout << ">" << seqID ++ << ": " << contigs[cp].contigID << ", " << contigs[cp].contigOffset << "| " << contigs[cp].contigID0 << ", " << contigs[cp].contigOffset0 << "|| " << contigs[cp].startID << ", " << contigs[cp].startOffset << "| " << contigs[cp].startID0 << ", " << contigs[cp].startOffset0 << "|| " << contigs[cp].endID << ", " << contigs[cp].endOffset << "| " << contigs[cp].endID0 << ", " << contigs[cp].endOffset0 << endl;
//			for(cpp = 0; cpp < contigs[cp].nucleotides.size(); cpp ++)
//			{
//				cout << contigs[cp].nucleotides[cpp];
//				if((cpp + 1) % 60 == 0 || cpp == contigs[cp].nucleotides.size() - 1)
//					cout << endl;
//			}
//		}	
//	}

}

void extendContigs(vector<vector<Base> > & genome, int coverage, int k, int chromosomeID)
{
	extdContigs1(genome, coverage, k, chromosomeID);
//	system("ps euf >> mem.txt");
//      genome.clear();
        extdContigs2(chromosomeID);
}

int overlap(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2) 
{
	if(x1 <= x2 && x2 <= y1 && y1 <= y2 && (int)y1 - (int)x2 > 0 || x2 <= x1 && x1 <= y2 && y2 <= y1 && (int)y2 - (int)x1 > 0 || x1 <= x2 && x2 <= y2 && y2 <= y1 && (int)y2 - (int)x2 > 0 || x2 <= x1 && x1 <= y1 && y1 <= y2 && (int)y1 - (int)x1 > 0)
		return 1;
	else
		return 0;
}

void scaffoldContigs(vector<vector<Base> > & genome, int chromosomeID)
{
	unsigned int seqID, cp, cp0, cpp, i, sp, spp;
	vector<vector<char> > scaffolds;
	vector<char> scaffold;
	int cont, covered;

	seqID = 0;
	for(cp = 0; cp < contigs.size(); cp ++)
	{
//		startIDBak = -1;
		if(contigs[cp].startID != -1 && contigs[cp].extended == 1)//added
		{
			scaffolds.push_back(scaffold);
			for(cpp = 0; cpp < contigs[cp].nucleotides.size(); cpp ++)
				scaffolds[scaffolds.size() - 1].push_back(contigs[cp].nucleotides[cpp]);
//			startIDBak = contigs[cp].startID;
			contigs[cp].startID = -1;
		}
		else
			continue;

		cont = 1;
		while(contigs[cp].startID0 == contigs[cp].endID0 && cont)
		{
			cont = 0;
			for(cp0 = cp + 1; cp0 < contigs.size(); cp0 ++)
			{
				if(cp0 != cp && contigs[cp].endID0 == contigs[cp0].startID && contigs[cp0].startID == contigs[cp0].endID && overlap(contigs[cp].startOffset0, contigs[cp].endOffset0, contigs[cp0].startOffset, contigs[cp0].endOffset) && contigs[cp0].extended == 1)//added
				{
					if(contigs[cp0].startOffset > contigs[cp].endOffset)
					{
						covered = 0;
						for(i = 0; i < contigs[cp0].startOffset - contigs[cp].endOffset - 1; i ++)
							if(genome[0][contigs[cp].endOffset + i + 1].kMer.size() > 0 || genome[0][contigs[cp].endOffset + i + 1].contiMer.size() > 0)
								covered ++;
						if(contigs[cp0].startOffset - contigs[cp].endOffset - 1 != 0 && (double) covered / (contigs[cp0].startOffset - contigs[cp].endOffset - 1) >= 0.5 || contigs[cp0].startOffset - contigs[cp].endOffset - 1 == 0)//0.8
							for(i = 0; i < contigs[cp0].startOffset - contigs[cp].endOffset - 1; i ++)
								scaffolds[scaffolds.size() - 1].push_back(genome[0][contigs[cp].endOffset + i + 1].nucleotide);
//								scaffolds[scaffolds.size() - 1].push_back('N');
						else
							goto conti;
					}
					for(cpp = 0; cpp < contigs[cp0].nucleotides.size(); cpp ++)
						scaffolds[scaffolds.size() - 1].push_back(contigs[cp0].nucleotides[cpp]);
//					startIDBak = contigs[cp0].startID;
					contigs[cp0].startID = -1;
					cp = cp0;
					cont = 1;
					break;
				}
conti:;
			}
		}	

	}

//	for(sp = 0; sp < scaffolds.size(); sp ++)
//	{
//		cout << ">" << chromosomeID << ": " << seqID ++ << endl;
//		for(spp = 0; spp < scaffolds[sp].size(); spp ++)
//		{
//			cout << scaffolds[sp][spp];
//			if((spp + 1) % 60 == 0 || spp == scaffolds[sp].size() - 1)
//				cout << endl;
//		}	
//	}

//output extended contigs
	ofstream out;
	string s = "tmp/_extended_contigs." + itoa(chromosomeID) + ".fa";
	out.open(s.c_str());
	for(sp = 0; sp < scaffolds.size(); sp ++)
	{
		out << ">" << seqID ++ << endl;
		for(spp = 0; spp < scaffolds[sp].size(); spp ++)
		{
			out << scaffolds[sp][spp];
			if((spp + 1) % 60 == 0 || spp == scaffolds[sp].size() - 1)
				out << endl;
		}
	}
}

int readLog()
{
        ifstream in;
        string buf;

        in.open("tmp/log.txt");
        if(in.is_open())
        {
                getline(in, buf);
                return(atoi(buf.c_str()));
        }
        else
        {
                cout << "CANNOT OPEN FILE!" << endl;
                exit(-1);
        }

}

void ref1()
{
	int numChromosomes, i, j, k, seqID;
	string s;
	vector<vector<char> > initContigs, extdContigs;
	vector<char> contig;
	string buf;
	unsigned int sourceID, targetID, targetStart, targetEnd, targetGap, sourceStart, sourceEnd, sourceGap, sourceSize, targetSize, fr;
	vector<int> initTags, extdTags;
	ifstream in, ex, ps;
	vector<Segment> seg;
	ofstream e, ini;
	vector<int> initNums;

	numChromosomes = readLog();

///////////////////////////////////////////////////////////////////for easy alignment/////////////////////////////////////////////////////////////////////////////////////

	ifstream tmpIn;
	ofstream tmpOut;
	vector<vector<char> > tmpContigs;
	for(i = 0; i < numChromosomes; i ++)
	{
		tmpIn.clear();
		tmpOut.clear();
		tmpContigs.clear();
		seqID = 0;

		s = "tmp/_initial_contigs." + itoa(i) + ".fa";
		tmpIn.open(s.c_str());
		s = "tmp/_short_initial_contigs." + itoa(i) + ".fa";
		tmpOut.open(s.c_str());
		if(tmpIn.is_open())
		{
			while(tmpIn.good())
			{
				getline(tmpIn, buf);
				if(buf[0] == 0) break;

				if(buf[0] == '>')
					tmpContigs.push_back(contig);
				else
					for(j = 0; j < buf.size(); j ++)
						tmpContigs[tmpContigs.size() - 1].push_back(buf[j]);
			}
		}
		else
		{
			cout << "CANNOT OPEN FILE!" << endl;
			exit(-1);
		}

		for(i = 0; i < tmpContigs.size(); i ++)
		{
			tmpOut << ">" << i << endl;
			if(tmpContigs[i].size() > SMALL_CHUNK)
			{
				for(j = 0; j < SMALL_CHUNK; j ++)
				{
					tmpOut << tmpContigs[i][j];
					if((j + 1) % 60 == 0 || j == SMALL_CHUNK - 1)//tmpContigs[i].size() - 1)
						tmpOut << endl;
				}
/*
				for(j = 20000; j < tmpContigs[i].size() - 20000; j ++)
				{
					tmpOut << "N";
					if((j + 1) % 60 == 0 || j == tmpContigs[i].size() - 1)
						tmpOut << endl;
				}
				for(j = tmpContigs[i].size() - 20000; j < tmpContigs[i].size(); j ++)
				{
					tmpOut << tmpContigs[i][j];
					if((j + 1) % 60 == 0 || j == tmpContigs[i].size() - 1)
						tmpOut << endl;
				}
*/
			}
			else
			{
				for(j = 0; j < tmpContigs[i].size(); j ++)
				{
					tmpOut << tmpContigs[i][j];
					if((j + 1) % 60 == 0 || j == tmpContigs[i].size() - 1)
						tmpOut << endl;
				}
			}
		}

		tmpIn.close();
		tmpOut.close();
	}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for(i = 0; i < numChromosomes; i ++)
	{
		s = "blat tmp/_extended_contigs." + itoa(i) + ".fa tmp/_short_initial_contigs." + itoa(i) + ".fa -noHead tmp/_short_initial_contigs_extended_contigs." + itoa(i) + ".psl  > blat_doc.txt";
		system(s.c_str());
	}
/*
	in.open("tmp/_contigs.fa");
	if(in.is_open())
	{
		while(in.good())
		{
			getline(in, buf);
			if(buf[0] == 0) break;

			if(buf[0] == '>')
			{
				initContigs.push_back(contig);
				initTags.push_back(0);
			}
			else
				for(i = 0; i < buf.size(); i++)
					initContigs[initContigs.size() - 1].push_back(buf[i]);
		}
	}
	else
	{
		cout << "CANNOT OPEN FILE!" << endl;
		return;
	}
*/
	seqID = 0;
	for(i = 0; i < numChromosomes; i ++)
	{
		s = "tmp/_extended_contigs." + itoa(i) + ".fa";
		ex.clear();
		ex.open(s.c_str());
		extdContigs.clear();
		extdTags.clear();
		if(ex.is_open())
		{
			while(ex.good())
			{
				getline(ex, buf);
				if(buf[0] == 0) break;

				if(buf[0] == '>')
				{
					extdContigs.push_back(contig);
					extdTags.push_back(0);
				}
				else
					for(j = 0; j < buf.size(); j ++)
						extdContigs[extdContigs.size() - 1].push_back(buf[j]);
			}
		}
		else
		{
			cout << "CANNOT OPEN FILE!" << endl;
			return;
		}

		s = "tmp/_initial_contigs." + itoa(i) + ".fa";
		in.clear();
	        in.open(s.c_str());
		initContigs.clear();
		initTags.clear();
	        if(in.is_open())
	        {
	                while(in.good())
	                {
	                        getline(in, buf);
	                        if(buf[0] == 0) break;

	                        if(buf[0] == '>')
	                        {
	                                initContigs.push_back(contig);
					initNums.push_back(atoi(buf.substr(1, buf.size()).c_str()));
	                                initTags.push_back(0);
	                        }
	                        else
	                                for(j = 0; j < buf.size(); j ++)
	                                        initContigs[initContigs.size() - 1].push_back(buf[j]);
	                }
	        }
	        else
	        {
	                cout << "CANNOT OPEN FILE!" << endl;
	                return;
	        }

		s = "tmp/_short_initial_contigs_extended_contigs." + itoa(i) + ".psl";
		ps.clear();
		ps.open(s.c_str());
		if(ps.is_open())
		{
			while(ps.good())
			{
				getline(ps, buf);
				if(buf[0] == 0) break;

				parseBLAT(buf, targetID, targetStart, targetEnd, targetGap, sourceID, sourceStart, sourceEnd, sourceGap, sourceSize, seg, fr, targetSize);
				if((double)(sourceEnd - sourceStart - sourceGap) / sourceSize >= 0.95 && (double)(targetEnd - targetStart - targetGap) / (double)(targetEnd - targetStart) >= 0.95 && targetSize > sourceSize + 100)
				{
					initTags[sourceID] = 1;
					extdTags[targetID] = 1;
				}
			}
		}
		else
		{
			cout << "CANNOT OPEN FILE!" << endl;
			return;
		}

		s = "tmp/_post_extended_contigs." + itoa(i) + ".fa";
		e.clear();
		e.open(s.c_str());
		if(e.is_open())
		{
			for(j = 0; j < extdTags.size(); j ++)
				if(extdTags[j] == 1)
				{
					e << ">" << seqID ++ << endl;
					for(k = 0; k < extdContigs[j].size(); k ++)
					{
						e << extdContigs[j][k];
						if(k == extdContigs[j].size() - 1 || (k + 1) % 60 == 0)
							e << endl;
					}
				}
		}
		else
		{
			cout << "CANNOT OPEN FILE!" << endl;
			return;
		}

		s = "tmp/_post_initial_contigs." + itoa(i) + ".fa";
		ini.clear();
		ini.open(s.c_str());
	        if(ini.is_open())
	        {
	                for(j = 0; j < initTags.size(); j ++)
	                        if(initTags[j] == 1)
	                        {
	                                ini << ">" << initNums[j] << endl;
	                                for(k = 0; k < initContigs[j].size(); k ++)
	                                {
	                                        ini << initContigs[j][k];
	                                        if(k == initContigs[j].size() - 1 || (k + 1) % 60 == 0)
	                                                ini << endl;
	                                }
	                        }
	        }
	        else
	        {
	                cout << "CANNOT OPEN FILE!" << endl;
	                return;
	        }

		ex.close();
		ps.close();
		e.close();
		ini.close();
	}
}

void ref2(ofstream & e, ofstream & r)
{
        int numChromosomes, i, j, k, seqID;
        string s;
        vector<vector<char> > initContigs, extdContigs;
        vector<char> contig;
        string buf;
        unsigned int sourceID, targetID, targetStart, targetEnd, targetGap, sourceStart, sourceEnd, sourceGap, sourceSize, targetSize, fr;
        vector<int> initTags, extdTags;
        ifstream in, ex, ps;
        vector<Segment> seg;

        numChromosomes = readLog();
        for(i = 0; i < numChromosomes; i ++)
        {
                s = "blat tmp/_post_extended_contigs." + itoa(i) + ".fa tmp/_post_initial_contigs." + itoa(i) + ".fa -noHead tmp/_initial_contigs_extended_contigs." + itoa(i) + ".psl > blat_doc.txt";
                system(s.c_str());
        }

        in.open("tmp/_contigs.fa");
        if(in.is_open())
        {
                while(in.good())
                {
                        getline(in, buf);
                        if(buf[0] == 0) break;

                        if(buf[0] == '>')
                        {
                                initContigs.push_back(contig);
                                initTags.push_back(0);
                        }
                        else
                                for(i = 0; i < buf.size(); i++)
                                        initContigs[initContigs.size() - 1].push_back(buf[i]);
                }
        }
        else
        {
                cout << "CANNOT OPEN FILE!" << endl;
                return;
        }

        seqID = 0;
        for(i = 0; i < numChromosomes; i ++)
        {
                s = "tmp/_post_extended_contigs." + itoa(i) + ".fa";
		ex.clear();
                ex.open(s.c_str());
                extdContigs.clear();
                extdTags.clear();
                if(ex.is_open())
                {
                        while(ex.good())
                        {
                                getline(ex, buf);
                                if(buf[0] == 0) break;

                                if(buf[0] == '>')
                                {
                                        extdContigs.push_back(contig);
                                        extdTags.push_back(0);
                                }
                                else
                                        for(j = 0; j < buf.size(); j ++)
                                                extdContigs[extdContigs.size() - 1].push_back(buf[j]);
                        }
                }
                else
                {
                        cout << "CANNOT OPEN FILE!" << endl;
                        return;
                }

                s = "tmp/_post_initial_contigs_post_extended_contigs." + itoa(i) + ".psl";
		ps.clear();
                ps.open(s.c_str());
                if(ps.is_open())
                {
                        while(ps.good())
                        {
                                getline(ps, buf);
                                if(buf[0] == 0) break;

                                parseBLAT(buf, targetID, targetStart, targetEnd, targetGap, sourceID, sourceStart, sourceEnd, sourceGap, sourceSize, seg, fr, targetSize);
                                if((double)(sourceEnd - sourceStart - sourceGap) / sourceSize >= 0.95 && (double)(targetEnd - targetStart - targetGap) / (double)(targetEnd - targetStart) >= 0.95 && targetSize > sourceSize + 100)
                                {
                                        initTags[sourceID] = 1;
                                        extdTags[targetID] = 1;
                                }
                        }
                }
                else
                {
                        cout << "CANNOT OPEN FILE!" << endl;
                        return;
                }

                if(e.is_open())
                {
                        for(j = 0; j < extdTags.size(); j ++)
                                if(extdTags[j] == 1)
                                {
                                        e << ">" << i << ": " << seqID ++ << endl;
                                        for(k = 0; k < extdContigs[j].size(); k ++)
                                        {
                                                e << extdContigs[j][k];
                                                if(k == extdContigs[j].size() - 1 || (k + 1) % 60 == 0)
                                                        e << endl;
                                        }
                                }
                }
                else
                {
                        cout << "CANNOT OPEN FILE!" << endl;
                        return;
                }
                ex.close();
                ps.close();
        }

        if(r.is_open())
        {
                for(i = 0; i < initTags.size(); i ++)
                        if(initTags[i] == 0)
                        {
                                r << ">" << i << endl;
                                for(j = 0; j < initContigs[i].size(); j ++)
                                {
                                        r << initContigs[i][j];
                                        if(j == initContigs[i].size() - 1 || (j + 1) % 60 == 0)
                                                r << endl;
                                }
                        }
        }
        else
        {
                cout << "CANNOT OPEN FILE!" << endl;
                return;
        }
}

//void refinement(ofstream & e, ofstream & r)
//{
//	ref1();
//	ref2(e, r);
//}

void refinement(ofstream & e, ofstream & r, int fastMap)
{
        int numChromosomes, i, j, k, seqID;
        string s, s0;
        vector<vector<char> > initContigs, extdContigs;
        vector<char> contig;
        string buf;
        unsigned int sourceID, targetID, targetStart, targetEnd, targetGap, sourceStart, sourceEnd, sourceGap, sourceSize, targetSize, fr, realSourceSize;
        vector<int> initTags, extdTags;
        ifstream in, ex, ps;
        vector<Segment> seg;
	vector<int> initNums;
	int ID, IDBak = -1;

#ifdef TEST
	ofstream ini, ext;
	ini.open("in.fa");
	ext.open("ex.fa");
#endif

        numChromosomes = readLog();

///////////////////////////////////////////////////////////////////for easy alignment/////////////////////////////////////////////////////////////////////////////////////

        ifstream tmpIn;
        ofstream tmpOut;
        vector<vector<char> > tmpContigs;
        for(i = 0; i < numChromosomes; i ++)
        {
                tmpIn.clear();
                tmpOut.clear();
                tmpContigs.clear();
		initNums.clear();
                seqID = 0;

                s = "tmp/_initial_contigs." + itoa(i) + ".fa";
                tmpIn.open(s.c_str());
                s = "tmp/_short_initial_contigs." + itoa(i) + ".fa";
                tmpOut.open(s.c_str());
                if(tmpIn.is_open())
                {
                        while(tmpIn.good())
                        {
                                getline(tmpIn, buf);
                                if(buf[0] == 0) break;

                                if(buf[0] == '>')
				{
                                        tmpContigs.push_back(contig);
					initNums.push_back(atoi(buf.substr(1, buf.size()).c_str()));
				}
                                else
                                        for(j = 0; j < buf.size(); j ++)
                                                tmpContigs[tmpContigs.size() - 1].push_back(buf[j]);
                        }
                }
                else
                {
                        cout << "CANNOT OPEN FILE!" << endl;
                        exit(-1);
                }

                for(j = 0; j < tmpContigs.size(); j ++)
                {
			
                        if(tmpContigs[j].size() > SMALL_CHUNK)
                        {
				tmpOut << ">" << initNums[j] << "." << tmpContigs[j].size()  << endl;// preserve id and size; parse for short contigs in a different way
                                for(k = 0; k < SMALL_CHUNK; k ++)
                                {
                                        tmpOut << tmpContigs[j][k];
                                        if((k + 1) % 60 == 0 || k == SMALL_CHUNK - 1)
                                                tmpOut << endl;
                                }
                        }
                        else
                        {
				tmpOut << ">" << initNums[j] << endl;
                                for(k = 0; k < tmpContigs[j].size(); k ++)
                                {
                                        tmpOut << tmpContigs[j][k];
                                        if((k + 1) % 60 == 0 || k == tmpContigs[j].size() - 1)
                                                tmpOut << endl;
                                }
                        }
                }

                tmpIn.close();
                tmpOut.close();
        }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(fastMap == 1) s0 = " -fastMap "; else s0 = " ";
        for(i = 0; i < numChromosomes; i ++)
        {
                s = "blat tmp/_extended_contigs." + itoa(i) + ".fa tmp/_short_initial_contigs." + itoa(i) + ".fa -noHead tmp/_short_initial_contigs_extended_contigs." + itoa(i) + ".psl" + s0 + ">> blat_doc.txt";
                system(s.c_str());
        }

        in.open("tmp/_contigs.fa");
        if(in.is_open())
        {
                while(in.good())
                {
                        getline(in, buf);
                        if(buf[0] == 0) break;

                        if(buf[0] == '>')
                        {
				for(i = 0; i < buf.size(); i ++)
					if(buf[i] == '.') break;
				ID = atoi(buf.substr(i + 1, buf.size()).c_str());
				if(ID != IDBak)
				{
                                	initContigs.push_back(contig);
                                	initTags.push_back(0);
					IDBak = ID;
				}
                        }
                        else
                                for(i = 0; i < buf.size(); i++)
                                        initContigs[initContigs.size() - 1].push_back(buf[i]);
                }
        }
        else
        {
                cout << "CANNOT OPEN FILE!" << endl;
                return;
        }

	int line;
	seqID = 0;
        for(i = 0; i < numChromosomes; i ++)
        {
                s = "tmp/_extended_contigs." + itoa(i) + ".fa";
                ex.open(s.c_str());
                extdContigs.clear();
                extdTags.clear();
                if(ex.is_open())
                {
                        while(ex.good())
                        {
                                getline(ex, buf);
                                if(buf[0] == 0) break;

                                if(buf[0] == '>')
                                {
                                        extdContigs.push_back(contig);
                                        extdTags.push_back(0);
                                }
                                else
                                        for(j = 0; j < buf.size(); j ++)
                                                extdContigs[extdContigs.size() - 1].push_back(buf[j]);
                        }
                }
                else
                {
                        cout << "CANNOT OPEN FILE!" << endl;
                        return;
                }

		line = 0;
                s = "tmp/_short_initial_contigs_extended_contigs." + itoa(i) + ".psl";
                ps.open(s.c_str());
                if(ps.is_open())
                {
                        while(ps.good())
                        {
                                getline(ps, buf);
				line ++;
                                if(buf[0] == 0) break;

                                realSourceSize = parseBLAT(buf, targetID, targetStart, targetEnd, targetGap, sourceID, sourceStart, sourceEnd, sourceGap, sourceSize, seg, fr, targetSize);
                                if((double)(sourceEnd - sourceStart - sourceGap) / sourceSize >= 0.8 && (double)(targetEnd - targetStart - targetGap) / (double)(targetEnd - targetStart) >= 0.8 && targetSize > realSourceSize + 100 && realSourceSize > targetSize / 100)
                                {
                                        initTags[sourceID] = 1;
                                        extdTags[targetID] = 1;
                                }
                        }
                }
                else
                {
                        cout << "CANNOT OPEN FILE!" << endl;
                        return;
                }

                if(e.is_open())
                {
                        for(j = 0; j < extdTags.size(); j ++)
                                if(extdTags[j] == 1)
                                {
                                        e << ">" << i << ": " << seqID << endl;
#ifdef TEST
					ext << ">" << i << ": " << seqID << endl;
#endif
					seqID ++;
                                        for(k = 0; k < extdContigs[j].size(); k ++)
                                        {
                                                e << extdContigs[j][k];
#ifdef TEST
						ext << extdContigs[j][k];
#endif
                                                if(k == extdContigs[j].size() - 1 || (k + 1) % 60 == 0)
						{
                                                        e << endl;
#ifdef TEST
							ext << endl;
#endif
						}
                                        }
                                }
                }
                else
                {
                        cout << "CANNOT OPEN FILE!" << endl;
                        return;
                }
                ex.close();
                ps.close();
        }

        if(r.is_open())
        {
                for(i = 0; i < initTags.size(); i ++)
                        if(initTags[i] == 0)
                        {
                                r << ">" << i << endl;
                                for(j = 0; j < initContigs[i].size(); j ++)
                                {
                                        r << initContigs[i][j];
                                        if(j == initContigs[i].size() - 1 || (j + 1) % 60 == 0)
                                                r << endl;
                                }
                        }
        }
        else
        {
                cout << "CANNOT OPEN FILE!" << endl;
                return;
        }

#ifdef TEST
        if(ini.is_open())
        {
                for(i = 0; i < initTags.size(); i ++)
                        if(initTags[i] == 1)
                        {
                                ini << ">" << i << endl;
                                for(j = 0; j < initContigs[i].size(); j ++)
                                {
                                        ini << initContigs[i][j];
                                        if(j == initContigs[i].size() - 1 || (j + 1) % 60 == 0)
                                                ini << endl;
                                }
                        }
        }
        else
        {
                cout << "CANNOT OPEN FILE!" << endl;
                return;
        }
#endif
}

int maxReadLength(ifstream & r)
{
	int max = 0, len = 0;
	string buf;

	if(r.is_open())
	{
		while(r.good())
		{
			getline(r, buf);
			if(buf[0] == 0) break;
			if(buf[0] == '>') 
			{
				if(len > max) 
					max = len; 
				len = 0; 
				continue;
			}
			len = len + buf.size();
		}
		if(len > max) max = len;
	}
	else
	{
		cout << "CANNOT OPEN FILE!" << endl;
		exit(-1);
	}

	return max;
}

void formalizeInput(ifstream & in, string file)
{
        string buf;
        int i, cp, cpp, total;
        unsigned long seqID = 0, count = 0, realSeqID = 0;
	vector<vector<char> > contigs;
	vector<char> contig;
        ofstream out;

	in.clear();
	in.seekg(0);
        out.open(file.c_str());

	if(file == "tmp/_contigs.fa")
	{
	        if(in.is_open())
	        {
	                while(in.good())
	                {
	                        getline(in, buf);
	                        if(buf[0] == 0)
	                                break;
	                        if(buf[0] == '>')
	                        {
					contigs.push_back(contig);
/*
					if(count % 60 != 0)
						out << endl;
					count = 0;
	                                in.getline(buf, 5000000);
	                                out << ">" << seqID ++ << endl;
	                                for(i = 0; i < in.gcount() - 1; i ++, count ++)
					{
	                                        out << buf[i];
						if((count + 1) % 60 == 0)
							out << endl;
					}
*/
	                        }
	                        else
	                        {
					for(i = 0; i < buf.size(); i ++)
						contigs[contigs.size() - 1].push_back(buf[i]);
/*
	                                for(i = 0; i < in.gcount() - 1; i ++, count ++)
					{
	                                        out << buf[i];
						if((count + 1) % 60 == 0)
	                                		out << endl;
					}
*/
	                        }
	                }

			for(cp = 0; cp < contigs.size(); cp ++)
			{
				if(contigs[cp].size() > 200)
				{
					if(contigs[cp].size() < LARGE_CHUNK)
					{
						out << ">" << seqID ++ << "." << realSeqID << endl;
						for(cpp = 0; cpp < contigs[cp].size(); cpp ++)
						{
							out << contigs[cp][cpp];
							if((cpp + 1) % 60 == 0 || cpp == contigs[cp].size() - 1)
								out << endl;
						}
					}
					else
					{
						out << ">" << seqID ++ << "." << realSeqID << endl;
						for(cpp = 0, total = 0; cpp < contigs[cp].size(); cpp ++)
						{
							out << contigs[cp][cpp];
							if((cpp + 1) % LARGE_CHUNK == 0 && cpp < contigs[cp].size() - 1 - 60)
							{
								total = total + LARGE_CHUNK;
								out << endl;
								out << ">" << seqID ++ << "." << realSeqID << endl;
								continue;
							}
							if((cpp + 1 - total) % 60 == 0 || cpp == contigs[cp].size() - 1)
								out << endl;
						}
					}
					realSeqID ++;
				}
			}
	        }
	        else
	        {
	                cout << "CANNOT OPEN FILE!" << endl;
	                exit(-1);
	        }
	}
	else
	{
		if(in.is_open())
		{
			while(in.good())
			{
				getline(in, buf);
				if(buf[0] == 0)
					break;
				if(buf[0] == '>')
					out << ">" << seqID ++ << endl;
				else
				{
					for(i = 0; i < buf.size(); i ++)
						out << buf[i];
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
}

int formalizeGenome(ifstream & in, int p)
{
        string buf;
        int i, chromosomeID, gp, cp, q;
        ofstream out;
	string s;
	vector<vector<char> > genome;
	vector<char> g;

        if(in.is_open())
        {
                while(in.good())
                {
                        getline(in, buf);
                        if(buf[0] == 0)
                                break;
                        if(buf[0] == '>')
                        {
				genome.push_back(g);
                        }
                        else
			{
                                for(i = 0; i < buf.size(); i ++)
                                        genome[genome.size() - 1].push_back(buf[i]);
			}
                }
        }
        else
        {
                cout << "CANNOT OPEN FILE!" << endl;
                exit(-1);
        }

	chromosomeID = 0;
	for(gp = 0; gp < genome.size(); gp ++)
	{
		s = "tmp/_genome." + itoa(chromosomeID) + ".fa";
		out.open(s.c_str());
		out << ">0" << endl;
		q = 1;
		for(cp = 0; cp < genome[gp].size(); cp ++)
		{
			out << genome[gp][cp];
			if((cp + 1) % 60 == 0 || cp == genome[gp].size() - 1 || ((cp + 1) % (genome[gp].size() / p) == 0 && q < p))
				out << endl;
			if(cp != genome[gp].size() - 1 && ((cp + 1) % (genome[gp].size() / p) == 0 && q < p))
			{
				out.close();
				chromosomeID ++;
				q ++;
				s = "tmp/_genome." + itoa(chromosomeID) + ".fa";
				out.open(s.c_str());
				out << ">0" << endl;
			}
		}
		out.close();
		chromosomeID ++;
	}
	genome.clear();
	out.close();
	return chromosomeID;
}

int formalizeInput(ifstream & in1, ifstream & in2, string file)
{
        string buf1, buf2;
        int i, count;
        unsigned long seqID = 0;
        ofstream out;

	in1.clear();
	in2.clear();
	in1.seekg(0);
	in2.seekg(0);
        out.open(file.c_str());
	count = 0;

        if(in1.is_open() && in2.is_open())
        {
                while(in1.good() && in2.good())
                {
                        getline(in1, buf1);
			getline(in2, buf2);
                        if(buf1[0] == 0 && buf2[0] == 0)
                                break;
			else if(buf1[0] == 0 && buf2[0] != 0 || buf1[0] != 0 && buf2[0] == 0)
			{
				cout << "INVALID INPUT FILE!" << endl;
				exit(-1);
			}

                        if(buf1[0] == '>' && buf2[0] == '>')
                        {
                                getline(in1, buf1);
				getline(in2, buf2);
                                out << ">" << seqID << endl;
                                for(i = 0; i < buf1.size(); i ++)
                                        out << buf1[i];
                                out << endl;
				out << ">" << seqID ++ << endl;
				for(i = 0; i < buf2.size(); i ++)
					out << buf2[i];
				out << endl;
				count ++;
                        }
                        else if(buf1[0] != '>' && buf2[0] != '>')
                        {
                                for(i = 0; i < buf1.size(); i ++)
                                        out << buf1[i];
                                out << endl;
                                for(i = 0; i < buf2.size(); i ++)
                                        out << buf2[i];
                                out << endl;
                        }
			else
			{
				cout << "INVALID INPUT FILE!" << endl;
				exit(-1);
			}
                }
        }
        else
        {
                cout << "CANNOT OPEN FILE!" << endl;
                exit(-1);
        }
	return count;
}

typedef struct insertStruct
{
	int distanceLow;
	int distanceHigh;
	int numChromosomes;
	int fastMap;
} Insert;

void * task0(void * arg)
{
	Insert ins;
	string command;
	stringstream distanceLowStr, distanceHighStr;
	int chromosomeID;

	ins = *(Insert *) arg;
	distanceLowStr << ins.distanceLow;
	distanceHighStr << ins.distanceHigh;
	ins.numChromosomes;
	ins.fastMap;

	for(chromosomeID = 0; chromosomeID < ins.numChromosomes; chromosomeID ++)
	{
		command = "bowtie2-build -f tmp/_genome." + itoa(chromosomeID) + ".fa tmp/_genome." + itoa(chromosomeID) + " > bowtie_doc.txt";
		system(command.c_str());
		command = "bowtie2 -f --no-mixed -k 5 -p 8 --local --mp 3,1 --rdg 2,1 --rfg 2,1 --score-min G,5,2 -I " + distanceLowStr.str() + " -X " + distanceHighStr.str() + " --no-discordant -x tmp/_genome." + itoa(chromosomeID) + " -1 tmp/_reads_1.fa -2 tmp/_reads_2.fa --reorder --very-sensitive-local > tmp/_reads_genome." + itoa(chromosomeID) + ".bowtie 2>> bowtie_doc.txt";// --very-sensitive-local
//                command = "bowtie2 -f --no-mixed -k 5 -p 8 --end-to-end --mp 3,1 --rdg 2,1 --rfg 2,1 --score-min L,-0.24,-0.24 -I " + distanceLowStr.str() + " -X " + distanceHighStr.str() + " --no-discordant -x tmp/_genome." + itoa(chromosomeID) + " -1 tmp/_reads_1.fa -2 tmp/_reads_2.fa --reorder > tmp/_reads_genome." + itoa(chromosomeID) + ".bowtie 2>> bowtie_doc.txt";
		system(command.c_str());
	}
}

void * task1(void * arg)
{
	Insert ins;
	int chromosomeID;
	string command, s;
	stringstream distanceLowStr, distanceHighStr;

	ins = *(Insert *) arg;
	distanceLowStr << ins.distanceLow;
	distanceHighStr << ins.distanceHigh;
	ins.numChromosomes;
	ins.fastMap;

	for(chromosomeID = 0; chromosomeID < ins.numChromosomes; chromosomeID ++)
	{
//		command = "lastdb -c -uMAM8 tmp/_genome." + itoa(chromosomeID) + ".db tmp/_genome." + itoa(chromosomeID) + ".fa > last_doc.txt";
//		system(command.c_str());
//		command = "lastal -e34 -m100 tmp/_genome." + itoa(chromosomeID) + ".db tmp/_contigs.fa | last-split > tmp/_contigs_genome." + itoa(chromosomeID) + ".maf 2>> last_doc.txt";
//		system(command.c_str());
//		command = "maf-convert.py psl tmp/_contigs_genome." + itoa(chromosomeID) + ".maf > tmp/_contigs_genome." + itoa(chromosomeID) + ".psl 2>> last_doc.txt";
//		system(command.c_str());
	}

	if(ins.fastMap == 1) s = " -fastMap "; else s = " ";
	for(chromosomeID = 0; chromosomeID < ins.numChromosomes; chromosomeID ++)
	{
		command = "blat tmp/_genome." + itoa(chromosomeID) + ".fa tmp/_contigs.fa -noHead tmp/_contigs_genome." + itoa(chromosomeID) + ".psl" + s + "> blat_doc.txt";// -minIdentity=50 -q=dnax -t=dnax
		system(command.c_str());
	}
}

void nonParallelMap(int distanceLow, int distanceHigh, int numChromosomes, int fastMap)
{
        string command, s;
        stringstream distanceLowStr, distanceHighStr;
	int chromosomeID;

	distanceLowStr << distanceLow;
	distanceHighStr << distanceHigh;
	if(fastMap == 1) s = " -fastmap "; else s = " ";
	for(chromosomeID = 0; chromosomeID < numChromosomes; chromosomeID ++)
	{
		command = "bowtie2-build -f tmp/_genome." + itoa(chromosomeID) + ".fa tmp/_genome." + itoa(chromosomeID) + " > bowtie_doc.txt";
		system(command.c_str());
		command = "bowtie2 -f --no-mixed -k 5 -p 8 --local --mp 3,1 --rdg 2,1 --rfg 2,1 --score-min G,5,2 -I " + distanceLowStr.str() + " -X " + distanceHighStr.str() + " --no-discordant -x tmp/_genome." + itoa(chromosomeID) + " -1 tmp/_reads_1.fa -2 tmp/_reads_2.fa --reorder --very-sensitive-local > tmp/_reads_genome." + itoa(chromosomeID) + ".bowtie 2>> bowtie_doc.txt";// --very-sensitive-local
//		command = "bowtie2 -f --no-mixed -k 5 -p 8 --end-to-end --mp 3,1 --rdg 2,1 --rfg 2,1 --score-min L,-0.24,-0.24 -I " + distanceLowStr.str() + " -X " + distanceHighStr.str() + " --no-discordant -x tmp/_genome." + itoa(chromosomeID) + " -1 tmp/_reads_1.fa -2 tmp/_reads_2.fa --reorder > tmp/_reads_genome." + itoa(chromosomeID) + ".bowtie 2>> bowtie_doc.txt";
		system(command.c_str());
		command = "blat tmp/_genome." + itoa(chromosomeID) + ".fa tmp/_contigs.fa -noHead tmp/_contigs_genome." + itoa(chromosomeID) + ".psl" + s + "> blat_doc.txt";// -minIdentity=50 -q=dnax -t=dnax
		system(command.c_str());
	}
}

void parallelMap(int distanceLow, int distanceHigh, int numChromosomes, int fastMap)
{
	pthread_t t0, t1;
	Insert ins;

	ins.distanceLow = distanceLow;
	ins.distanceHigh = distanceHigh;
	ins.numChromosomes = numChromosomes;
	ins.fastMap = fastMap;
	if(pthread_create(&t0, NULL, task0, &ins) != 0) {nonParallelMap(distanceLow, distanceHigh, numChromosomes, fastMap); return;}
	if(pthread_create(&t1, NULL, task1, &ins) != 0) {nonParallelMap(distanceLow, distanceHigh, numChromosomes, fastMap); return;}
	
	if(pthread_join(t0, NULL) != 0) {perror("Thread join faied"); exit(EXIT_FAILURE);}
	if(pthread_join(t1, NULL) != 0) {perror("Thread join faied"); exit(EXIT_FAILURE);}
}

void writeLog(int n)
{
	ofstream out;

	out.open("tmp/log.txt");
	if(out.is_open())
		out << n << endl;
	else
	{
		cout << "CANNOT OPEN FILE!" << endl;
		exit(-1);
	}
}

void checkRatio(int numChromosomes)
{
	int i, rp, chromosomeID, aligned;
	unsigned targetID1, targetStart1, targetEnd1, targetGap1, sourceID1, sourceStart1, sourceEnd1, sourceGap1, sourceSize1, fr1, targetID2, targetStart2, targetEnd2, targetGap2, sourceID2, sourceStart2, sourceEnd2, sourceGap2, sourceSize2, fr2;
        vector<Segment> segs1, segs2;
	ifstream r, ra;
	string s, buf;
	vector<int> reads;
	double ratio;

	r.open("tmp/_reads_1.fa");
	if(r.is_open())
	{
		while(r.good())
		{
			getline(r, buf);
			if(buf[0] == '>')
				reads.push_back(0);
		}
	}
	else
	{
		cout << "CANNOT OPEN FILE!" << endl;
		exit(-1);
	}

	for(chromosomeID = 0; chromosomeID < numChromosomes; chromosomeID ++)
	{
		s = "tmp/_reads_genome." + itoa(chromosomeID) + ".bowtie";
		ra.open(s.c_str());
		if(ra.is_open())
		{
			while(ra.good())
			{
				getline(ra, buf);
				if(buf[0] == 0) break;
				if(buf[0] == '@') continue;
				parseBOWTIE(buf, targetID1, targetStart1, targetEnd1, targetGap1, sourceID1, sourceStart1, sourceEnd1, sourceGap1, sourceSize1, segs1, fr1);

				getline(ra, buf);
				if(buf[0] == 0)
				{
					cout << "BROKEN BOWTIE FILE" << endl;
					exit(-1);
				}
				parseBOWTIE(buf, targetID2, targetStart2, targetEnd2, targetGap2, sourceID2, sourceStart2, sourceEnd2, sourceGap2, sourceSize2, segs2, fr2);
				if(targetID1 != -1 && targetID2 != -1 && (double) (sourceEnd1 - sourceStart1 - sourceGap1) / /*(sourceEnd1 - sourceStart1)*/ sourceSize1 >= THRESHOLD && (double) (targetEnd1 - targetStart1 - targetGap1) / (targetEnd1 - targetStart1) >= THRESHOLD && (double) (sourceEnd2 - sourceStart2 - sourceGap2) / /*(sourceEnd2 - sourceStart2)*/ sourceSize2 >= THRESHOLD && (double) (targetEnd2 - targetStart2 - targetGap2) / (targetEnd2 - targetStart2) >= THRESHOLD)
					reads[sourceID1] = 1;
			}
		}
		else
		{
			cout << "CANNOT OPEN FILE!" << endl;
			exit(-1);
		}
	}

	for(rp = 0, aligned = 0; rp < reads.size(); rp ++)
		if(reads[rp] == 1)
			aligned ++;
	ratio = aligned == 0 ? 0 : (double) aligned / reads.size();
	cout << " - " << ratio * 100 << "% reads aligned ";
	if(ratio < 0.5 && ratio >= 0.25)
		cout << "(warning: ratio a little low; may not guarantee good results)" << endl;
	else if(ratio < 0.25)
		cout << "(warning: ratio too low; hard to guarantee good results)" << endl;
	else
		cout << endl;
}

void print()
{
	cout << "AlignGraph --read1 reads_1.fa --read2 reads_2.fa --contig contigs.fa --genome genome.fa --distanceLow distanceLow --distanceHigh distancehigh --extendedContig extendedContigs.fa --remainingContig remainingContigs.fa [--kMer k --insertVariation insertVariation --covereage coverage --noAlignment --part p --checkRatio]" << endl;
	cout << "Inputs:" << endl;
	cout << "--read1 is the the first pair of PE DNA reads in fasta format" << endl;
	cout << "--read2 is the the second pair of PE DNA reads in fasta format" << endl;
	cout << "--contig is the initial contigs in fasta format" << endl;
	cout << "--genome is the reference genome in fasta format" << endl;
	cout << "--distanceLow is the lower bound of alignment distance between the first and second pairs of PE DNA reads" << endl;
	cout << "--distanceHigh is the upper bound of alignment distance between the first and second pairs of PE DNA reads" << endl;
	cout << "Outputs:" << endl;
	cout << "--extendedContig is the extended contig file in fasta format" << endl;
	cout << "--remainingContig is the not extended initial contig file in fasta format" << endl;
	cout << "Options:" << endl;
	cout << "--kMer is the k-mer size (default: 5)" << endl;
	cout << "--insertVariation is the small variation of insert length (default: 50)" << endl;
	cout << "--coverage is the minimum coverage to keep a path in de Bruijn graph (default: 20)" << endl;
	cout << "--noAlignment skips the initial time-consuming alignment step, if all the alignment files have been provided in tmp directory (default: none)" << endl;
	cout << "--part is the number of parts a chromosome is divided into when it is loaded to reduce memory requirement (default: 1)" << endl;
	cout << "--fastMap makes BLAT alignment faster but may lower LEAF's performance. Useful for large genomes (default: none)" << endl;
	cout << "--checkRatio checks read alignment ratio to the reference beforehand and warns if the ratio is too low; may take a little more time (default: none)" << endl;
}

int main(int argc, char * argv[])
{
	ifstream r1, r2, c, g;
	ofstream e, r;
	string s, sCheck;
	stringstream ssCheck;
	int i, tagRead1 = 0, tagRead2 = 0, tagContig = 0, tagGenome = 0, tagExtendedContig = 0, tagKMer = 0, tagDistanceLow = 0, tagDistanceHigh = 0, tagNoAlignment = 1, k = 5, distanceLow = 0, distanceHigh = MAX, chromosomeID, numChromosomes, coverage = 20, tagCoverage = 0, mrl, mrl1, mrl2, tagInsertVariation = 0, insertVariation = 50, tagRemainingContig = 0, part = 1, tagPart = 0, tagFastMap = 0, numReads, tagCheckRatio = 0;
	time_t start, end, startAlign, endAlign;

//	e.open("extended_contigs.fa");
//	r.open("remaining_contigs.fa");
//	refinement(e, r, 0);
//	return 1;

	cout << "AlignGraph: algorithm for secondary de novo genome assembly guided by closely related references" << endl;
	cout << "By Ergude Bao, CS Department, UC-Riverside. All Rights Reserved" << endl << endl;

	start = time(NULL);
	for(i = 1; i < argc; i ++)
	{
		s = argv[i];
                if(s == "--read1")
                {
                        if(tagRead1 == 1 || i == argc - 1)
                        {
                                print();
                                return 0;
                        }
                        r1.open(argv[++ i]);
                        if(!r1.is_open())
                        {
                                cout << "CANNOT OPEN FILE!" << endl;
                                print();
                                return 0;
                        }
			tagRead1 = 1;
                }
                else if(s == "--read2")
                {
                        if(tagRead2 == 1 || i == argc - 1)
                        {
                                print();
                                return 0;
                        }
                        r2.open(argv[++ i]);
                        if(!r2.is_open())
                        {
                                cout << "CANNOT OPEN FILE!" << endl;
                                print();
                                return 0;
                        }
                        tagRead2 = 1;
                }
                else if(s == "--contig")
                {
                        if(tagContig == 1 || i == argc - 1)
                        {
                                print();
                                return 0;
                        }
                        c.open(argv[++ i]);
                        if(!c.is_open())
                        {
                                cout << "CANNOT OPEN FILE!" << endl;
                                print();
                                return 0;
                        }
                        tagContig = 1;
                }
                else if(s == "--genome")
                {
                        if(tagGenome == 1 || i == argc - 1)
                        {
                                print();
                                return 0;
                        }
                        g.open(argv[++ i]);
                        if(!g.is_open())
                        {
                                cout << "CANNOT OPEN FILE!" << endl;
                                print();
                                return 0;
                        }
                        tagGenome = 1;
                }
                else if(s == "--distanceLow")
                {
                        if(tagDistanceLow == 1 || i == argc - 1)
                        {
                                print();
                                return 0;
                        }
                        distanceLow = atoi(argv[++ i]);
                        ssCheck.str("");
                        ssCheck << distanceLow;
                        sCheck = ssCheck.str();
                        if(sCheck != argv[i])
                        {
                                print();
                                return 0;
                        }
                        tagDistanceLow = 1;
                }
                else if(s == "--distanceHigh")
                {
                        if(tagDistanceHigh == 1 || i == argc - 1)
                        {
                                print();
                                return 0;
                        }
                        distanceHigh = atoi(argv[++ i]);
                        ssCheck.str("");
                        ssCheck << distanceHigh;
                        sCheck = ssCheck.str();
                        if(sCheck != argv[i])
                        {
                                print();
                                return 0;
                        }
                        tagDistanceHigh = 1;
                }
                else if(s == "--extendedContig")
                {
                        if(tagExtendedContig == 1 || i == argc - 1)
                        {
                                print();
                                return 0;
                        }
			e.open(argv[++ i]);
                        if(!e.is_open())
                        {
                                cout << "CANNOT OPEN FILE!" << endl;
                                print();
                                return 0;
                        }
			tagExtendedContig = 1;
                }
                else if(s == "--remainingContig")
                {
                        if(tagRemainingContig == 1 || i == argc - 1)
                        {
                                print();
                                return 0;
                        }
                        r.open(argv[++ i]);
                        if(!r.is_open())
                        {
                                cout << "CANNOT OPEN FILE!" << endl;
                                print();
                                return 0;
                        }
                        tagRemainingContig = 1;
                }
                else if(s == "--kMer")
                {
                        if(tagKMer == 1 || i == argc - 1)
                        {
                                print();
                                return 0;
                        }
			k = atoi(argv[++ i]);
			ssCheck.str("");
			ssCheck << k;
			sCheck = ssCheck.str();
			if(sCheck != argv[i])
			{
				print();
				return 0;
			}
			tagKMer = 1;
                }
                else if(s == "--insertVariation")
                {
                        if(tagInsertVariation == 1 || i == argc - 1)
                        {
                                print();
                                return 0;
                        }
                        insertVariation = atoi(argv[++ i]);
                        ssCheck.str("");
                        ssCheck << insertVariation;
                        sCheck = ssCheck.str();
                        if(sCheck != argv[i])
                        {
                                print();
                                return 0;
                        }
                        tagInsertVariation = 1;
                }
                else if(s == "--coverage")
                {
                        if(tagCoverage == 1 || i == argc - 1)
                        {
                                print();
                                return 0;
                        }
                        coverage = atoi(argv[++ i]);
                        ssCheck.str("");
                        ssCheck << coverage;
                        sCheck = ssCheck.str();
                        if(sCheck != argv[i])
                        {
                                print();
                                return 0;
                        }
                        tagCoverage = 1;
                }
                else if(s == "--noAlignment")
                {
                        if(tagNoAlignment == 0)
                        {
                                print();
                                return 0;
                        }
                        tagNoAlignment = 0;
                }
		else if(s == "--part")
		{
			if(tagPart == 1 || i == argc - 1)
			{
				print();
				return 0;
			}
			part = atoi(argv[++ i]);
			ssCheck.str("");
			ssCheck << part;
			sCheck = ssCheck.str();
			if(sCheck != argv[i])
			{
				print();
				return 0;
			}
			tagPart = 1;
		}
		else if(s == "--fastMap")
		{
			if(tagFastMap == 1)
			{
				print();
				return 0;
			}
			tagFastMap = 1;
		}
		else if(s == "--checkRatio")
		{
			if(tagCheckRatio == 1)
			{
				print();
				return 0;
			}
			tagCheckRatio = 1;
		}
		else
		{
			print();
			return 0;
		}
	}

	if(tagRead1 == 0 || tagRead2 == 0 || tagContig == 0 || tagGenome == 0 || tagExtendedContig == 0 || tagRemainingContig == 0 || k <= 0 || tagDistanceLow == 0 || tagDistanceHigh == 0 || distanceLow > distanceHigh || distanceLow < 0 || insertVariation < 0 || part < 1 || part > 10)// || k > mrl)
	{
		print();
		return 0;
	}
	mrl1 = maxReadLength(r1);
	mrl2 = maxReadLength(r2);
	mrl = mrl1 > mrl2 ? mrl1 : mrl2;
	if(k > mrl)
	{
		print();
		return 0;
	}

	system("test -d \"tmp\"; t=$?; if [ $t -eq 1 ]; then mkdir tmp; fi");

	if(tagNoAlignment == 1)
	{
		formalizeInput(r1, r2, "tmp/_reads.fa");
		formalizeInput(r1, "tmp/_reads_1.fa");
		formalizeInput(r2, "tmp/_reads_2.fa");
		formalizeInput(c, "tmp/_contigs.fa");
		numChromosomes = formalizeGenome(g, part);
		startAlign = time(NULL);
		parallelMap(distanceLow, distanceHigh, numChromosomes, tagFastMap);
		endAlign = time(NULL);
		writeLog(numChromosomes);
		cout << "(0) Alignment finished" << endl;
	}
	else
	{
		numChromosomes = readLog();
		startAlign = endAlign = time(NULL);
		cout << "(0) Alignment skipped" << endl;
	}

	if(tagCheckRatio == 1)
		checkRatio(numChromosomes);

	for(chromosomeID = 0; chromosomeID < numChromosomes; chromosomeID ++)
	{
		cout << endl << "CHROMOSOME " << chromosomeID << ": " << endl;
		loadGenome(genome, chromosomeID);
		cout << "(1) chromosome loaded" << endl;
		loadContigAlignment(genome, chromosomeID);
		cout << "(2) contig alignment loaded" << endl;
		loadReadAlignment(genome, k, insertVariation, chromosomeID, mrl);
		cout << "(3) read alignment loaded" << endl;
		extendContigs(genome, coverage, k, chromosomeID);
		cout << "(4) contigs extended" << endl;
		scaffoldContigs(genome, chromosomeID);
		cout << "(5) contigs scaffolded" << endl;
		system("ps euf >> mem.txt");
		contigs.clear();
		genome.clear();
		sourceIDBak = -1;
	}

	refinement(e, r, tagFastMap);

	end = time(NULL);
	cout << endl << "FINISHED for " <<  end - start << " seconds (" << endAlign - startAlign << " seconds for alignment) :-)" << endl;
}

