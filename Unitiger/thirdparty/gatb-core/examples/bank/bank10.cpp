//! [snippet1]
// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>
#include <sstream>

// We use the required packages
using namespace std;

static const char* STR_BANKS_NB       = "-split";
static const char* STR_MAX_INPUT_SIZE = "-max-size";

/********************************************************************************/
/*                         Bank split                                           */
/*                                                                              */
/* This snippet shows how to split a bank into smaller banks and how to create  */
/* an album bank (ie a list of URL of banks). Such an album bank could be used  */
/* as bank input by other tools.                                                */
/* Note: all the generated files are put in a directory created by the snippet. */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    /** We create a command line parser. */
    OptionsParser parser ("BankSplitter");
    parser.push_back (new OptionOneParam (STR_URI_INPUT,      "bank reference",            true));
    parser.push_back (new OptionOneParam (STR_MAX_INPUT_SIZE, "average db size per split", true));
    parser.push_back (new OptionOneParam (STR_BANKS_NB,       "number max of sub banks",   false, "0"));
    parser.push_back (new OptionOneParam (STR_URI_OUTPUT_DIR, "output directory",          false, "."));

    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        /** We parse the user options. */
        IProperties* options = parser.parse (argc, argv);

        /** Shortcuts. */
        u_int64_t maxDbSize = options->getInt(STR_MAX_INPUT_SIZE);

        // We declare an input Bank
        IBank* inputBank = BankRegistery::singleton().createBank (options->getStr(STR_URI_INPUT));
        LOCAL (inputBank);

        // We get the basename of the input bank.
        string inputBasename = System::file().getBaseName (options->getStr(STR_URI_INPUT));

        int nbBanks = options->getInt (STR_BANKS_NB);

        /** We set the name of the output directory. */
        stringstream ss;  ss << inputBasename << "_S" << maxDbSize << "_N" << nbBanks;
        string outputDirName = ss.str();

        /** We create the output directory. */
        string outputDir = options->getStr(STR_URI_OUTPUT_DIR) + "/" + outputDirName;
        System::file().mkdir (outputDir, S_IRWXU);

        // We create the album bank.
        BankAlbum album (outputDir + "/album.txt");

        // We create a sequence iterator for the bank
        Iterator<Sequence>* itInput = inputBank->iterator();

        /** We get estimations about the bank. */
        u_int64_t number, totalSize, maxSize;
        inputBank->estimate (number, totalSize, maxSize);

        u_int64_t estimationNbSeqToIterate = nbBanks <= 0 ? number  : (number*maxDbSize*nbBanks) / totalSize;

        // We create an iterator over the input bank and encapsulate it with progress notification.
        SubjectIterator<Sequence> itSeq (itInput, 1000);
        itSeq.addObserver (new ProgressTimer (estimationNbSeqToIterate, "split"));

        // We loop over sequences to get the exact number of sequences.
          int64_t nbBanksOutput = -1;
        u_int64_t nbSequences   =  0;
        u_int64_t dbSize        = ~0;

        IBank* currentBank = 0;

        for (itSeq.first(); !itSeq.isDone(); itSeq.next())
        {
            if (dbSize > maxDbSize)
            {
                if (currentBank != 0)  { currentBank->flush(); }

                nbBanksOutput ++;
                if (nbBanks > 0 && nbBanksOutput >= nbBanks)  { break; }

                /** We build the uri of the current bank. */
                stringstream ss;  ss << inputBasename << "_" << nbBanksOutput;

                /** We create a new bank and put it in the album. */
                currentBank = album.addBank (outputDir, ss.str());

                /** We reinit the db size counter. */
                dbSize = 0;
            }

            dbSize += itSeq->getDataSize();

            /** We insert the sequence into the current output bank. */
            currentBank->insert (*itSeq);
        }

        if (currentBank != 0)  { currentBank->flush(); }
    }
    catch (OptionFailure& e)
    {
        e.getParser().displayErrors (stdout);
        e.getParser().displayHelp   (stdout);
        return EXIT_FAILURE;
    }
    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }
}
//! [snippet1]
