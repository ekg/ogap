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


short qualityChar2ShortInt(char c) {
    return static_cast<short>(c) - 33;
}

void countMismatchesAndGaps(
    BamAlignment& alignment,
    vector<CigarOp>& cigarData,
    string referenceSequence,
    int& mismatches,
    int& gaps,
    int& gapslen,
    int& softclips,
    int& mismatchQsum,
    int& softclipQsum
    ) {

    int sp = 0;
    int rp = 0;
    for (vector<CigarOp>::const_iterator c = cigarData.begin();
	 c != cigarData.end(); ++c) {
        int l = c->Length;
        char t = c->Type;
        if (t == 'M') { // match or mismatch
            for (int i = 0; i < l; ++i) {
                if (alignment.QueryBases.at(rp) != referenceSequence.at(sp)) {
                    ++mismatches;
		    mismatchQsum += qualityChar2ShortInt(alignment.Qualities.at(rp));
		}
                ++sp;
                ++rp;
            }
        } else if (t == 'D') { // deletion
            ++gaps;
	    gapslen += l;
            sp += l;  // update reference sequence position
        } else if (t == 'I') { // insertion
	    ++gaps;
	    gapslen += l;
	    rp += l;  // update read position
	} else if (t == 'S') { // soft clip, clipped sequence present in the read not matching the reference
	    softclips += l;
	    for (int i = 0; i < l; ++i) {
		softclipQsum += qualityChar2ShortInt(alignment.Qualities.at(rp));
		++rp;
	    }
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
	 << "    -x --mismatch-rate N       Trigger realignment if the mismatch rate is greater" << endl
	 << "                               than N (default 0.03)" << endl
	 << "    -c --mismatch-count N      Trigger realignment if the mismatch count is greater" << endl
	 << "                               than or equal to N (default -1, unbounded)" << endl
	 << "    -q --soft-clip-count N     Trigger realignment if there are more than N softclips" << endl
	 << "                               (default -1, unbounded)" << endl
	 << "    -C --mismatch-qsum N       Trigger realignment if the sum of quality scores of" << endl
	 << "                               mismatched bases is >= N" << endl
	 << "    -Q --soft-clip-qsum N      Trigger realignment if the sum of quality scores of" << endl
	 << "                               soft clipped bases is >= N" << endl
	 << "    -S --soft-clip-limit N     Only accept realignments if they have <= N soft clips" << endl
	 << "    -i --max-gap-increase N    Only allow the introduction of up to N gaps" << endl
	 << "                               soft-clipped bases is >= N (default 3)" << endl
	 << "    -M --min-mapping-quality N Only realign reasd with MQ > N (default 0)" << endl
	 << "    -Z --zero-mapping-quality  Set mapping quality to 0 when the read would be realigned" << endl
	 << "    -m --match-score N         The match score (default 10.0)" << endl
	 << "    -n --mismatch-score N      The mismatch score (default -9.0)" << endl
	 << "    -g --gap-open-penalty N    The gap open penalty (default 15.0)" << endl
	 << "    -e --gap-extend-penalty N  The gap extend penalty (default 6.66)" << endl
	 << "    -z --entropy-gap-open     use entropy scaling for the gap open penalty" << endl
	 << "    -R --repeat-gap-extend N  penalize non-repeat-unit gaps in repeat sequence" << endl
	 << "                               (default 15.0, set to 0 to disable)" << endl
	 << "    -s --suppress-output       Don't output BAM on stdout" << endl
	 << "    -A --accept-all            Accept all realignments, regardless of quality" << endl;   
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

    float mismatchRateAfterLimit = 0.1f;

    int softclipCountThreshold = 0;

    int mismatchQualitySumThreshold = -1;
    int softclipQualitySumThreshold = -1;

    bool acceptAllRealignments = false;

    float acceptGapsToSoftclips = 5.0f;
    float acceptGapsToMismatches = 1.0f;

    int maxGapIncrease = 3;

    int minMappingQuality = 0;

    int softclipLimit = -1;

    bool zeroMq = false;
    
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
	    {"mismatch-rate", required_argument, 0, 'x'},
	    {"mismatch-count", required_argument, 0, 't'},
	    {"soft-clip-count", required_argument, 0, 'q'},
	    {"mismatch-qsum", required_argument, 0, 'C'},
	    {"soft-clip-qsum", required_argument, 0, 'Q'},
            {"match-score",  required_argument, 0, 'm'},
            {"mismatch-score",  required_argument, 0, 'n'},
            {"gap-open-penalty",  required_argument, 0, 'g'},
            {"gap-extend-penalty",  required_argument, 0, 'e'},
	    {"entropy-gap-open", no_argument, 0, 'z'},
	    {"repeat-gap-extend", required_argument, 0, 'R'},
            {"suppress-output", no_argument, 0, 's'},
	    {"accept-all", no_argument, 0, 'A'},
	    {"max-gap-increase", required_argument, 0, 'i'},
	    {"min-mapping-quality", required_argument, 0, 'M'},
	    {"soft-clip-limit", required_argument, 0, 'S'},
	    {"zero-mapping-quality", no_argument, 0, 'Z'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hszZAdf:m:n:g:e:w:c:x:R:q:C:Q:i:M:S:",
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

            case 'q':
		softclipCountThreshold = atoi(optarg);
		break;

            case 'S':
		softclipLimit = atoi(optarg);
		break;

	    case 'Z':
		zeroMq = true;
		break;

	    case 'C':
		mismatchQualitySumThreshold = atoi(optarg);
		break;

	    case 'M':
		minMappingQuality = atoi(optarg);
		break;

	    case 'A':
		acceptAllRealignments = true;
		break;

            case 'Q':
		softclipQualitySumThreshold = atoi(optarg);
		break;

	    case 'i':
		maxGapIncrease = atoi(optarg);
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
        if (alignment.IsMapped() && alignment.MapQuality >= minMappingQuality) {

	    string ref;

	    try {

		/*
	    int softclipBegin = 0;
	    int softclipEnd = 0;
	    if (alignment.CigarData.front().Type == Constants::BAM_CIGAR_SOFTCLIP_CHAR)
		softclipBegin = alignment.CigarData.front().Length;
	    if (alignment.CigarData.back().Type == Constants::BAM_CIGAR_SOFTCLIP_CHAR)
		softclipEnd = alignment.CigarData.back().Length;

            int endpos = alignment.GetEndPosition() + softclipEnd;
            int length = endpos - (alignment.Position - softclipBegin);  // 0-based half-open interval
		*/

            int endpos = alignment.GetEndPosition();
            int length = endpos - alignment.Position;  // 0-based half-open interval

            // do we meet criteria for realignment?

	    // get the overlapping reference sequnce to determine mismatches
	    ref = reference.getSubSequence(referenceIDToName[alignment.RefID],
					   max(0, alignment.Position - flanking_window),
					   length + 2 * flanking_window);

	    if (debug) cerr << ref << endl;
	    if (debug) cerr << alignment.QueryBases << endl;

	    int mismatchesBefore = 0;
	    int gapsBefore = 0;
	    int gapslenBefore = 0;
	    int softclipsBefore = 0;
	    int mismatchQsumBefore = 0;
	    int softclipQsumBefore = 0;
	    countMismatchesAndGaps(alignment,
				   alignment.CigarData,
				   ref.substr(flanking_window, length),
				   mismatchesBefore,
				   gapsBefore,
				   gapslenBefore,
				   softclipsBefore,
				   mismatchQsumBefore,
				   softclipQsumBefore);


	    // if we meet the criteria, attempt to realign

            if (gapsBefore > 0
		|| (softclipCountThreshold > -1 && softclipsBefore > softclipCountThreshold)
		|| (mismatchCountThreshold > -1 && mismatchesBefore >= mismatchCountThreshold)
		|| (softclipQualitySumThreshold > -1 && softclipQsumBefore >= softclipQualitySumThreshold)
		|| (mismatchQualitySumThreshold > -1 && mismatchQsumBefore >= mismatchQualitySumThreshold)
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
		int gapslenAfter = 0;
		int softclipsAfter = 0;
		int mismatchQsumAfter = 0;
		int softclipQsumAfter = 0;
		countMismatchesAndGaps(alignment,
				       cigarData,
				       ref.substr(referencePos, realignmentLength),
				       mismatchesAfter,
				       gapsAfter,
				       gapslenAfter,
				       softclipsAfter,
				       mismatchQsumAfter,
				       softclipQsumAfter);


		if (debug) {
		    cerr << cigar << endl;
		    cerr << alignment.QueryBases << endl;
		    cerr << ref.substr(referencePos, realignmentLength) << endl;
		    cerr << "before: " << mismatchesBefore << " mismatches, " << gapsBefore << " gaps, and " << softclipsBefore << " softclips" << endl;
		    cerr << "after:  " << mismatchesAfter << " mismatches, " << gapsAfter << " gaps, and " << softclipsAfter << " softclips" << endl;
		    cerr << "alignment delta: " << alignmentStartDelta << " length delta " << alignmentLengthDelta << " end delta " << alignmentEndDelta << endl;
		}

		     /*
		     && (!(softclipsBefore > softclipsAfter && gapsBefore < gapsAfter)
			 || ((softclipsBefore > softclipsAfter && gapsBefore < gapsAfter) &&
			 (float) (softclipsBefore - softclipsAfter)
			     / (float) (gapsBefore - gapsAfter) > acceptGapsToSoftclips)
		     && (float) variancesAfter
			 / (float) alignment.QueryBases.size() < mismatchRateAfterLimit
		     */



		// reduce mismatches or maintain them,
		// reduce the sum of variances in the alignment
		//int variancesBefore = mismatchesBefore + gapsBefore + softclipsBefore;
		//int variancesAfter = mismatchesAfter + gapsAfter + softclipsAfter;
		//if (debug) cerr << mismatchQsumAfter + softclipQsumAfter << " ? < " << mismatchQsumBefore + softclipQsumBefore << endl;
		if (acceptAllRealignments ||
		    (mismatchQsumAfter <= mismatchQsumBefore
		     && softclipQsumAfter <= softclipQsumBefore
		     //(mismatchQsumAfter + softclipQsumAfter <= mismatchQsumBefore + softclipQsumBefore
		     && (softclipLimit == -1 || (softclipLimit >= 0 && softclipsAfter <= softclipLimit))
		     && ((gapsBefore == gapsAfter && gapslenAfter <= gapslenBefore)
			 || gapsAfter - gapsBefore <= maxGapIncrease)
		     //&& variancesAfter <= variancesBefore
		     //&& mismatchesAfter + gapsAfter <= mismatchesBefore + gapsBefore
		     //&& (alignmentStartDelta <= flanking_window  // don't accept -pos alignments
		     //    && alignmentEndDelta <= flanking_window
		     //    && alignment.Position + alignmentStartDelta >= 0))) {
			)) {

		    if (debug) cerr << "adjusting ... " << endl;
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
		    if (zeroMq) {
			alignment.MapQuality = 0;
		    } else {
			alignment.CigarData = cigarData;
			alignment.Position += alignmentStartDelta;
		    }
                }

		if (debug) cerr << endl;
            }

	} catch (...) {
		cerr << "exception when realigning " << alignment.Name << " at position " << referenceIDToName[alignment.RefID] << ":" << alignment.Position << " " << ref << " " << alignment.QueryBases << endl;
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
