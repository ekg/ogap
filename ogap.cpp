#include <iostream>
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <signal.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <map>
#include <vector>

#include "Fasta.h"
#include "BamAlignment.h"
#include "BamReader.h"
#include "BamWriter.h"

#include "SmithWatermanGotoh.h"

using namespace std;
using namespace BamTools;


pair<int, int> countMismatchesAndGaps(BamAlignment& alignment, string referenceSequence) {

    pair<int, int> mgps = make_pair(0, 0);
    int sp = 0;
    int rp = 0;
    for (vector<CigarOp>::const_iterator c = alignment.CigarData.begin();
        c != alignment.CigarData.end(); ++c) {
        int l = c->Length;
        char t = c->Type;
        if (t == 'M') { // match or mismatch
            for (int i = 0; i < l; ++i) {
                if (alignment.QueryBases.at(rp) != referenceSequence.at(sp))
                    ++mgps.first;
                ++sp;
                ++rp;
            }
        } else if (t == 'D') { // deletion
            ++mgps.second;
            sp += l;  // update reference sequence position
        } else if (t == 'I') { // insertion
            ++mgps.second;
            rp += l;  // update read position
        } else if (t == 'S') { // soft clip, clipped sequence present in the read not matching the reference
            rp += l;
        } else if (t == 'H') { // hard clip on the read, clipped sequence is not present in the read
        } else if (t == 'N') { // skipped region in the reference not present in read, aka splice
            sp += l;
        }
    }

    return mgps;

}


void printUsage(char** argv) {
    cerr << "usage: [BAM data stream] | " << argv[0] << " [options]" << endl
         << endl
         << "Realigns alignments meeting specified criteria (number of gaps, mismatches) using" << endl
         << "Smith-Waterman parameters optimized to open gaps and eliminate mismatches and" << endl
         << "writes the stream of alignments as BAM on stdout." << endl
         << endl
         << "arguments:" << endl
         << "    -f --fasta-reference FILE  FASTA reference file to use for realignment (required)" << endl
         << "    -d --debug                 Print debugging information about realignment process" << endl
         << "    -w --flanking-window BP    The number of bases on the left and right to attempt to" << endl
         << "                               align to (default 300bp)." << endl
         << "    -c --max-position-delta    Maximum number of bases the start or end of the alignment" << endl
         << "                               may change when realigning (default 200bp)" << endl
         //<< "    -z --max-mismatch-delta    Maximum change in number of mismatches" << endl
         << "    -m --match-score           The match score (default 10.0)" << endl
         << "    -n --mismatch-score        The mismatch score (default -20.0)" << endl
         << "    -g --gap-open-penalty      The gap open penalty (default 15.0)" << endl
         << "    -e --gap-extend-penalty    The gap extend penalty (default 0.0)" << endl
         << "    -s --suppress-output       Don't output BAM on stdout" << endl;
}

int main(int argc, char** argv) {

    int c;

    FastaReference reference;
    bool has_ref = false;
    bool suppress_output = false;
    bool debug = false;
    int flanking_window = 200;

    float matchScore = 10.0f;
    float mismatchScore = -20.0f;
    float gapOpenPenalty = 15.0f;
    float gapExtendPenalty = 0.0f;

    int maxPositionDelta = 200;
    
    if (argc < 2) {
        printUsage(argv);
        exit(1);
    }

    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"debug", no_argument, 0, 'd'},
            {"fasta-reference", required_argument, 0, 'f'},
            {"flanking-window", required_argument, 0, 'w'},
            {"max-position-delta", required_argument, 0, 'c'},
            {"match-score",  required_argument, 0, 'm'},
            {"mismatch-score",  required_argument, 0, 'n'},
            {"gap-open-penalty",  required_argument, 0, 'g'},
            {"gap-extend-penalty",  required_argument, 0, 'e'},
            {"suppress-output", no_argument, 0, 's'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hsdf:m:n:g:e:w:c:",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;
 
        switch (c) {

            case 'f':
                reference.open(optarg); // will exit on open failure
                has_ref = true;
                break;
     
            case 'd':
                debug = true;
                break;

            case 'w':
                flanking_window = atoi(optarg);
                break;

            case 'c':
                maxPositionDelta = atoi(optarg);
                break;

            case 'm':
                matchScore = atof(optarg);
                break;
 
            case 'n':
                mismatchScore = atof(optarg);
                break;
 
            case 'g':
                gapOpenPenalty = atof(optarg);
                break;
 
            case 'e':
                gapExtendPenalty = atof(optarg);
                break;

            case 's':
                suppress_output = true;
                break;

            case 'h':
                printUsage(argv);
                exit(0);
                break;
              
            case '?':
                printUsage(argv);
                exit(1);
                break;
     
              default:
                abort();
                break;
        }
    }

    if (!has_ref) {
        cerr << "no FASTA reference provided, cannot realign" << endl;
        exit(1);
    }

    BamReader reader;
    if (!reader.Open("stdin")) {
        cerr << "could not open stdin for reading" << endl;
        exit(1);
    }

    BamWriter writer;
    if (!suppress_output && !writer.Open("stdout", reader.GetHeaderText(), reader.GetReferenceData(), true)) {
        cerr << "could not open stdout for writing" << endl;
        exit(1);
    }

    // store the names of all the reference sequences in the BAM file
    map<int, string> referenceIDToName;
    vector<RefData> referenceSequences = reader.GetReferenceData();
    int i = 0;
    for (RefVector::iterator r = referenceSequences.begin(); r != referenceSequences.end(); ++r) {
        referenceIDToName[i] = r->RefName;
        ++i;
    }

    BamAlignment alignment;

    while (reader.GetNextAlignment(alignment)) {
        
        // skip unmapped alignments, as they cannot be left-realigned without CIGAR data
        if (alignment.IsMapped()) {

            int endpos = alignment.GetEndPosition();
            int length = endpos - alignment.Position + 1;

            // do we meet criteria for realignment?
            // ...
            // is there an indel?
            bool hasindel = false;

            stringstream cigar_before, cigar_after;
            for (vector<CigarOp>::const_iterator c = alignment.CigarData.begin();
                c != alignment.CigarData.end(); ++c) {
                unsigned int l = c->Length;
                char t = c->Type;
                cigar_before << l << t;
                if (t == 'D' || t == 'I') {
                    hasindel = true;
                    break;
                }
            }

            // if not, continue
            if (hasindel) {

                const string ref = reference.getSubSequence(referenceIDToName[alignment.RefID], alignment.Position - flanking_window, length + 2 * flanking_window);

                pair<int, int> mgps_before = countMismatchesAndGaps(alignment, ref.substr(flanking_window, length));

                const unsigned int referenceLen = ref.size();
                const unsigned int queryLen     = alignment.QueryBases.size();

                unsigned int referencePos;
                string cigar;

                // realign

                CSmithWatermanGotoh sw(matchScore, mismatchScore, gapOpenPenalty, gapExtendPenalty);
                sw.Align(referencePos, cigar, ref.c_str(), referenceLen, alignment.QueryBases.c_str(), queryLen);

                // parse the cigar

                vector<CigarOp> cigarData;

                string slen;
                int len = 0;
                int realignmentLength = 0;

                for (string::iterator c = cigar.begin(); c != cigar.end(); ++c) {
                    switch (*c) {
                        case 'I':
                            len = atoi(slen.c_str());
                            slen.clear();
                            cigarData.push_back(CigarOp(*c, len));
                            break;
                        case 'D':
                            len = atoi(slen.c_str());
                            realignmentLength += len;
                            slen.clear();
                            cigarData.push_back(CigarOp(*c, len));
                            break;
                        case 'M':
                            len = atoi(slen.c_str());
                            realignmentLength += len;
                            slen.clear();
                            cigarData.push_back(CigarOp(*c, len));
                            break;
                        case 'S':
                            len = atoi(slen.c_str());
                            //referencePos += len;
                            //realignmentLength += len;
                            slen.clear();
                            cigarData.push_back(CigarOp(*c, len));
                            break;
                        case 'N':
                            len = atoi(slen.c_str());
                            slen.clear();
                            cigarData.push_back(CigarOp(*c, len));
                            break;
                        case 'H':
                            len = atoi(slen.c_str());
                            //referencePos += len;  ?
                            slen.clear();
                            cigarData.push_back(CigarOp(*c, len));
                            break;
                        default:
                            len = 0;
                            slen += *c;
                            break;
                    }
                }

                int alignmentStartDelta = referencePos - flanking_window;
                int alignmentLengthDelta = realignmentLength - length;
                int alignmentEndDelta = alignmentStartDelta + alignmentLengthDelta;

                // if we have a delta, but it's within bounds, modify the alignment
                pair<int, int> mgps_after = countMismatchesAndGaps(alignment, ref.substr(flanking_window, length));

                if (((mgps_after.first < mgps_before.first
                        && mgps_after.second <= mgps_before.second)
                        ||
                    (mgps_after.first <= mgps_before.first
                        && mgps_after.second < mgps_before.second))
                    && (alignmentStartDelta <= maxPositionDelta
                        && alignmentEndDelta <= maxPositionDelta)) {
                    alignment.CigarData = cigarData;
                    alignment.Position += alignmentStartDelta;
                }
            }
        }

        // write every alignment unless we are suppressing output
        if (!suppress_output)
            writer.SaveAlignment(alignment);

    }

    reader.Close();
    if (!suppress_output)
        writer.Close();

    return 0;

}
