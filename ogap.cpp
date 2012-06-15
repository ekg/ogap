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

#include "fastahack/Fasta.h"

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAlignment.h"

#include "smithwaterman/SmithWatermanGotoh.h"

using namespace std;
using namespace BamTools;


void countMismatchesAndGaps(BamAlignment& alignment, vector<CigarOp>& cigarData, string referenceSequence, int& mismatches, int& gaps, int& softclips) {

    int sp = 0;
    int rp = 0;
    for (vector<CigarOp>::const_iterator c = cigarData.begin();
        c != cigarData.end(); ++c) {
        int l = c->Length;
        char t = c->Type;
        if (t == 'M') { // match or mismatch
            for (int i = 0; i < l; ++i) {
                if (alignment.QueryBases.at(rp) != referenceSequence.at(sp))
                    ++mismatches;
                ++sp;
                ++rp;
            }
        } else if (t == 'D') { // deletion
            ++gaps;
            sp += l;  // update reference sequence position
        } else if (t == 'I') { // insertion
	    ++gaps;
            rp += l;  // update read position
        } else if (t == 'S') { // soft clip, clipped sequence present in the read not matching the reference
	    softclips += l;
            rp += l;
        } else if (t == 'H') { // hard clip on the read, clipped sequence is not present in the read
        } else if (t == 'N') { // skipped region in the reference not present in read, aka splice
            sp += l;
        }
    }

}


void printUsage(char** argv) {
    cerr << "usage: [BAM data stream] | " << argv[0] << " [options]" << endl
         << endl
         << "Realigns alignments meeting specified criteria (number of gaps, mismatches) using" << endl
         << "a repeat-aware Smith-Waterman alignment algorithm optimized to open gaps in low-entropy" << endl
	 << "sequence.  New alignments are kept if the alignment reduces the overall number of" << endl
	 << "mismatches, gaps, and soft-clipped bases relative to the initial alignment while not" << endl
	 << "increasing the number of mismatches.  Realigned alignments are written as BAM on stdout." << endl
         << endl
         << "arguments:" << endl
         << "    -f --fasta-reference FILE  FASTA reference file to use for realignment (required)" << endl
         << "    -d --debug                 Print debugging information about realignment process" << endl
         << "    -w --flanking-window N     The number of bases on the left and right to attempt to" << endl
         << "                               align to (default 50bp).  This limits the maximum detectable" << endl
	 << "                               indel length." << endl
	 << "    -x --mismatch-rate-threshold N   Trigger realignment if the mismatch rate is greater" << endl
	 << "                                     than N (default 0.03)" << endl
	 << "    -c --mismatch-count-threshold N  Trigger realignment if the mismatch count is greater" << endl
	 << "                                     than or equal to N (default -1, unbounded)" << endl
         << "    -m --match-score N         The match score (default 10.0)" << endl
         << "    -n --mismatch-score N      The mismatch score (default -9.0)" << endl
         << "    -g --gap-open-penalty N    The gap open penalty (default 15.0)" << endl
         << "    -e --gap-extend-penalty N  The gap extend penalty (default 6.66)" << endl
	 << "    -z, --entropy-gap-open     use entropy scaling for the gap open penalty" << endl
	 << "    -R, --repeat-gap-extend N  penalize non-repeat-unit gaps in repeat sequence" << endl
	 << "                               (default 15.0, set to 0 to disable)" << endl
         << "    -s --suppress-output       Don't output BAM on stdout" << endl;
}

int main(int argc, char** argv) {

    int c;

    FastaReference reference;
    bool has_ref = false;
    bool suppress_output = false;
    bool debug = false;
    int flanking_window = 50;

    float matchScore = 10.0f;
    float mismatchScore = -9.0f;
    float gapOpenPenalty = 15.0f;
    float gapExtendPenalty = 6.66f;

    bool useEntropy = false;
    bool useRepeatGapExtendPenalty = true;
    float repeatGapExtendPenalty = 15.0f;

    float mismatchRateThreshold = 0.03f;
    int mismatchCountThreshold = -1;
    
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
	    {"mismatch-rate-threshold", required_argument, 0, 'x'},
	    {"mismatch-count-threshold", required_argument, 0, 't'},
            {"match-score",  required_argument, 0, 'm'},
            {"mismatch-score",  required_argument, 0, 'n'},
            {"gap-open-penalty",  required_argument, 0, 'g'},
            {"gap-extend-penalty",  required_argument, 0, 'e'},
	    {"entropy-gap-open", no_argument, 0, 'z'},
	    {"repeat-gap-extend", no_argument, 0, 'R'},
            {"suppress-output", no_argument, 0, 's'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hszdf:m:n:g:e:w:c:x:R:",
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

	    case 'R':
		repeatGapExtendPenalty = atof(optarg);
		if (repeatGapExtendPenalty == 0) {
		    useRepeatGapExtendPenalty = false;
		}
		break;

	    case 'z':
	        useEntropy = true;
		break;

	    case 'x':
		mismatchRateThreshold = atof(optarg);
		break;

	    case 'c':
		mismatchCountThreshold = atoi(optarg);
		break;

            case 'w':
                flanking_window = atoi(optarg);
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
    if (!suppress_output && !writer.Open("stdout", reader.GetHeaderText(), reader.GetReferenceData())) {
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
    map<long unsigned int, vector<BamAlignment> > alignmentSortQueue;

    while (reader.GetNextAlignment(alignment)) {

	long unsigned int initialAlignmentPosition = alignment.Position;
        
	// todo, align unmapped reads
        if (alignment.IsMapped()) {

            int endpos = alignment.GetEndPosition();
            int length = endpos - alignment.Position + 1;

            // do we meet criteria for realignment?

	    // get the overlapping reference sequnce to determine mismatches
	    const string ref = reference.getSubSequence(referenceIDToName[alignment.RefID],
							max(0, alignment.Position - flanking_window),
							length + 2 * flanking_window);

	    int mismatchesBefore = 0;
	    int gapsBefore = 0;
	    int softclipsBefore = 0;
	    countMismatchesAndGaps(alignment,
				   alignment.CigarData,
				   ref.substr(flanking_window, length),
				   mismatchesBefore,
				   gapsBefore,
				   softclipsBefore);

	    // if we meet the criteria, attempt to realign

            if (gapsBefore > 0
		|| softclipsBefore > 0
		|| (mismatchCountThreshold > -1 && mismatchesBefore >= mismatchCountThreshold)
		|| (float) mismatchesBefore / (float) alignment.QueryBases.size() > mismatchRateThreshold
		) {

                unsigned int referencePos;
                string cigar;

                // realign

                CSmithWatermanGotoh sw(matchScore, mismatchScore, gapOpenPenalty, gapExtendPenalty);
		if (useEntropy) sw.EnableEntropyGapPenalty(1);
		if (useRepeatGapExtendPenalty) sw.EnableRepeatGapExtensionPenalty(repeatGapExtendPenalty);
                sw.Align(referencePos, cigar, ref, alignment.QueryBases);

                // parse the cigar

                vector<CigarOp> cigarData;

                string slen;
                int len = 0;
                int realignmentLength = 0;
		int deletedBases = 0;

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
			    deletedBases += len;
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
                            realignmentLength += len;
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

		// todo, deal with weird alignment artifact: deletion
		// at the front of a read without any matching,
		// anchoring, sequence

		int alignmentStartDelta;
		if (alignment.Position - flanking_window < 0) {
		    alignmentStartDelta = (int) referencePos - alignment.Position;
		} else {
		    alignmentStartDelta = (int) referencePos - flanking_window;
		}
                int alignmentLengthDelta = realignmentLength - length;
                int alignmentEndDelta = alignmentStartDelta + alignmentLengthDelta;

                // if we have a delta, but it's within bounds, modify the alignment

		int mismatchesAfter = 0;
		int gapsAfter = 0;
		int softclipsAfter = 0;
		countMismatchesAndGaps(alignment,
				       cigarData,
				       ref.substr(referencePos, realignmentLength),
				       mismatchesAfter,
				       gapsAfter,
				       softclipsAfter);
		/*
		cerr << cigar << endl;
		cerr << alignment.QueryBases << endl;
		cerr << ref.substr(referencePos, realignmentLength) << endl;
		cerr << "before: " << mismatchesBefore << " mismatches, " << gapsBefore << " gaps, and " << softclipsBefore << " softclips" << endl;
		cerr << "after:  " << mismatchesAfter << " mismatches, " << gapsAfter << " gaps, and " << softclipsAfter << " softclips" << endl;
		cerr << "alignment delta: " << alignmentStartDelta << " length delta " << alignmentLengthDelta << " end delta " << alignmentEndDelta << endl;
		*/

		// reduce mismatches or maintain them,
		// reduce the sum of variances in the alignment
		if (mismatchesAfter <= mismatchesBefore
		    && mismatchesAfter + gapsAfter + softclipsAfter
		    <= mismatchesBefore + gapsBefore + softclipsBefore
                    && (alignmentStartDelta <= flanking_window  // don't accept -pos alignments
                        && alignmentEndDelta <= flanking_window
			&& alignment.Position + alignmentStartDelta >= 0)) {
		    //cerr << "adjusting ... " << endl;
		    /*
		    if (abs(alignmentStartDelta) > deletedBases) { // weird... why?
			cerr << "odd move of " << alignmentStartDelta << " for " << alignment.Name << endl;
			for (vector<CigarOp>::iterator o = cigarData.begin(); o != cigarData.end(); ++o)
			    cerr << o->Length << o->Type;
			cerr << endl;
			cerr << alignment.Position << endl;
			cerr << ref.substr(referencePos, realignmentLength) << endl;
			cerr << alignment.AlignedBases << endl << endl;
	            }
		    */
		    alignment.CigarData = cigarData;
		    alignment.Position += alignmentStartDelta;
                }

		//cerr << endl;
            }
        }

        // write every alignment unless we are suppressing output
        if (!suppress_output) {
	    alignmentSortQueue[alignment.Position].push_back(alignment);
	    // ensure correct order if alignments move
	    if (initialAlignmentPosition > (unsigned int) flanking_window) {
		long unsigned int maxOutputPos = initialAlignmentPosition - flanking_window;
		map<long unsigned int, vector<BamAlignment> >::iterator p = alignmentSortQueue.begin();
		for ( ; p != alignmentSortQueue.end(); ++p) {
		    if (p->first > maxOutputPos) {
			break; // no more to do
		    } else {
			for (vector<BamAlignment>::iterator a = p->second.begin(); a != p->second.end(); ++a)
			    writer.SaveAlignment(*a);
		    }
		}
		if (p != alignmentSortQueue.begin())
		    alignmentSortQueue.erase(alignmentSortQueue.begin(), p);
	    }
	}

    }

    if (!suppress_output) {
	for (map<long unsigned int, vector<BamAlignment> >::iterator p = alignmentSortQueue.begin();
	     p != alignmentSortQueue.end(); ++p) {
	    for (vector<BamAlignment>::iterator a = p->second.begin(); a != p->second.end(); ++a)
		writer.SaveAlignment(*a);
	}
	alignmentSortQueue.clear();
    }


    reader.Close();
    if (!suppress_output)
        writer.Close();

    return 0;

}
