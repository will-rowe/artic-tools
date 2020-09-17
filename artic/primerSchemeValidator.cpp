#include <fstream>

#include "log.hpp"
#include "primerScheme.hpp"

// ValidateScheme will load and validate a specified primer scheme, log stats and return the primer scheme object.
// It will optionally write primer insert coordinates and primer sequences to files.
artic::PrimerScheme artic::ValidateScheme(SchemeArgs& args)
{
    // check there is a scheme loaded or try loading one
    LOG_TRACE("reading scheme")
    if (args.schemeFile.size() == 0)
        LOG_ERROR("no primer scheme file was provided");
    auto ps = artic::PrimerScheme(args.schemeFile, args.schemeVersion);

    // get primer sequences if requested
    if (args.primerSeqsFile.size() != 0)
    {
        LOG_TRACE("collecting primer sequences")
        if (args.refSeqFile.size() == 0)
            LOG_ERROR("no reference sequence provided, can't output primer sequences");
        faidx_t* fai = fai_load(args.refSeqFile.c_str());
        std::ofstream fh;
        fh.open(args.primerSeqsFile);
        for (auto amplicon : ps.GetExpAmplicons())
        {
            auto f = (amplicon.GetForwardPrimer()->GetNumAlts()) ? amplicon.GetForwardPrimer()->GetID() + std::string("_alts_merged") : amplicon.GetForwardPrimer()->GetID();
            auto r = (amplicon.GetReversePrimer()->GetNumAlts()) ? amplicon.GetReversePrimer()->GetID() + std::string("_alts_merged") : amplicon.GetReversePrimer()->GetID();
            fh << ">" << f << std::endl;
            fh << amplicon.GetForwardPrimer()->GetSeq(fai, ps.GetReferenceName()) << std::endl;
            fh << ">" << r << std::endl;
            fh << amplicon.GetReversePrimer()->GetSeq(fai, ps.GetReferenceName()) << std::endl;
        }
        fh.close();
        if (fai)
            fai_destroy(fai);
        LOG_TRACE("\twritten to file: {}", args.primerSeqsFile);
    }

    // get primer insert coordinates if requested
    if (args.insertsFile.size() != 0)
    {
        LOG_TRACE("collecting primer insert coordinates")
        std::ofstream fh;
        fh.open(args.insertsFile);
        int counter = 1;
        for (auto amplicon : ps.GetExpAmplicons())
        {
            auto poolID = amplicon.GetPrimerPoolID();
            fh << ps.GetReferenceName() << "\t" << amplicon.GetForwardPrimer()->GetEnd() << "\t" << amplicon.GetReversePrimer()->GetStart() << "\t" << counter << "\t" << ps.GetPrimerPool(poolID) << "\t+" << std::endl;
            counter++;
        }
        fh.close();
        LOG_TRACE("\twritten to file: {}", args.insertsFile);
    }

    // print the stats
    LOG_TRACE("collecting scheme stats");
    LOG_TRACE("\tprimer scheme version:\t{}", ps.GetVersion());
    LOG_TRACE("\treference sequence:\t{}", ps.GetReferenceName());
    LOG_TRACE("\tnumber of pools:\t{}", ps.GetPrimerPools().size());
    LOG_TRACE("\tnumber of primers:\t{} (includes {} alts)", ps.GetNumPrimers(), ps.GetNumAlts());
    LOG_TRACE("\tnumber of amplicons:\t{}", ps.GetNumAmplicons());
    LOG_TRACE("\tmean amplicon size:\t{}", ps.GetMeanAmpliconSpan());
    LOG_TRACE("\tscheme ref. span:\t{}-{}", ps.GetRefStart(), ps.GetRefEnd());
    float proportion = (float)ps.GetNumOverlaps() / (float)(ps.GetRefEnd() - ps.GetRefStart());
    LOG_TRACE("\tscheme overlaps:\t{}%", proportion * 100);
    return ps;
}